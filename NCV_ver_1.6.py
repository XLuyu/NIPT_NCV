# coding: utf-8
import os
import time
import argparse
import pickle
from datetime import datetime
import pysam
import pp
import pandas as pd
from tqdm.autonotebook import tqdm
import genome_browser

def tqdm_pp_jobs(jobs, desc="pp jobs"):
    result = [tqdm_pp_jobs]*len(jobs)
    with tqdm(total=len(jobs), desc=desc) as pbar:
        while tqdm_pp_jobs in result: 
            for i in range(len(result)):
                if result[i]==tqdm_pp_jobs and jobs[i].finished:
                    result[i] = jobs[i]()
                    pbar.update(1)
            time.sleep(1)
            pbar.update(0)
    return result

def get_readcount_table(filepath):
    bamfile = pysam.AlignmentFile(filepath, "rb")
    chroms = ["chr" + str(i) for i in range(1,23)] + ["chrX", "chrY"]
    chrlen = {i['SN']:i['LN'] for i in bamfile.header['SQ'] if i['SN'] in chroms}
    count_array = {}
    for chr in chroms:
        array = [0]*(chrlen[chr]//50000+1)
        for read in bamfile.fetch(chr):
            pos = read.reference_end-1 if read.is_reverse else read.reference_start
            array[pos//50000] += 1
        count_array[chr] = array
    return count_array

######## Train ########
def optimize(coi, training):  # the index of chromesome of interest(1-based)
    import itertools
    chroms = ["chr" + str(i) for i in range(1, 23)]
    if coi in chroms: chroms.remove(coi)
    optimized = [float("inf")]
    for i in range(1, len(chroms) + 1):
        for comb in itertools.combinations(chroms, i):
            comb = list(comb)
            ratios = training[coi] / training[comb].sum(axis=1)
            mean = ratios.mean()
            std = ratios.std()
            CV = std / mean
            if CV < optimized[0]:
                optimized = [CV, coi, comb, mean, std]
    return optimized[1:]

def NCVtrain(args):
    ### input
    M_files = [os.path.abspath(line.strip()) for line in open(args.male) if line.strip() != ""]
    F_files = [os.path.abspath(line.strip()) for line in open(args.female) if line.strip() != ""]
    job_server = pp.Server(min(6, args.cpus), ppservers=("*:3456",))
    M_counts = tqdm_pp_jobs([job_server.submit(get_readcount_table, (f,), modules=("pysam",)) for f in M_files], desc="Male train samples")
    F_counts = tqdm_pp_jobs([job_server.submit(get_readcount_table, (f,), modules=("pysam",)) for f in F_files], desc="Female train samples")
    chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    df = pd.DataFrame({chr:[sum(counts[chr]) for counts in M_counts+F_counts] for chr in chroms}, index=M_files+F_files, dtype="float")
    df['Gender'] = ["M"]*len(M_counts) + ["F"]*len(F_counts)
    ### train reference chromosome set
    inputs = [("chr"+str(i), df) for i in range(1, 23)]
    inputs += [("chrX", df[df.Gender=="F"]), ("chrY", df[df.Gender=="F"])]
    job_server = pp.Server(args.cpus, ppservers=("*:3456",))
    refset = tqdm_pp_jobs([job_server.submit(optimize, i, modules=("numpy", "pandas")) for i in inputs], desc="Enumerating reference set")
    ### calculate mean/std for each bin
    bin_mean_std = {}
    pickle.dump((refset, bin_mean_std), open(args.model, "w"))
    for chr, comb, _, _ in tqdm(refset, desc="Bin training"):
        bin_num = [len(counts[chr]) for counts in M_counts+F_counts]
        if min(bin_num)!=max(bin_num): raise Exception("BAMs from different ref genome")
        bin_mean_std[chr] = []
        for i in range(bin_num[0]):
            bin_count = [counts[chr][i] for counts in (F_counts if chr in ['chrX', 'chrY'] else M_counts+F_counts)]
            ref_count = df[df.Gender=="F"] if chr in ['chrX', 'chrY'] else df
            bin_ratio = bin_count/ref_count[comb].sum(axis=1)
            bin_mean_std[chr] += [(bin_ratio.mean(), bin_ratio.std())]
    pickle.dump((refset, bin_mean_std), open(args.model, "w"))
            

######## Test ########
def getXY(filename): # get read counts on chrX and chrY
    def check(read): return not read.is_unmapped and read.mapping_quality >= 60
    bam = pysam.AlignmentFile(filename,"rb")
    tot = float(bam.count(read_callback=check))
    return bam.count('chrY', read_callback=check)/tot#, bam.count('chrX', read_callback=check)/tot

def FF(X_count, Y_count):
    fx,fy,my,mx = 0.0509404113275, 2.32850602014e-05, 0.00212467247788, 0.0288937505057
    FF_x = (X_count-fx)/(mx-fx)
    FF_y = (Y_count-fy)/(my-fy)
    return FF_x, FF_y

def buildNCVTableTests(counts, refset):
    zscore = pd.DataFrame({},index=counts.index)
    for numer, denom, mean, std in refset:
        zscore[numer] = (counts[numer] / counts[denom].sum(axis=1) - mean) / std
    return zscore

def NCVtest(args):
    ### input
    refset, bin_mean_std = pickle.load(open(args.model, "rb"))
    filenames = [os.path.abspath(line.strip()) for line in open(args.test) if line.strip() != ""]
    job_server = pp.Server(min(6, args.cpus), ppservers=("*:3456",))
    test = tqdm_pp_jobs([job_server.submit(get_readcount_table, (f,), modules=("pysam",)) for f in filenames], desc="Test samples")
    chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    df = pd.DataFrame({chr:[sum(counts[chr]) for counts in test] for chr in chroms}, index=filenames, dtype="float")
    ### test for chromosome
    result_trisomy = buildNCVTableTests(df, refset)
    result_trisomy = result_trisomy[chroms]
    fx, fy, my, mx = 0.0509404113275, 2.32850602014e-05, 0.00212467247788, 0.0288937505057
    result_trisomy['FF'] = tqdm_pp_jobs([job_server.submit(getXY, (f,), modules=("pysam",)) for f in filenames], desc="FF test")
    result_trisomy['FF'] = (result_trisomy['FF']-fy)/(my-fy)
    os.mkdir(args.output)
    result_trisomy.to_html(open(args.output + "/table.html", "w"))
    result_trisomy.to_csv(open(args.output + "/table.csv", "w"))
    ### test for micro-insertion/deletion
    gb = genome_browser.GenomeBroswer()
    for f, counts in tqdm(zip(filenames, test), desc="Chromosome plot per sample"):
        fname = os.path.basename(f)[:-4]
        os.mkdir(args.output+"/"+fname)
        for chr, comb, _, _ in refset:
            x, y = [], []
            for i in range(len(bin_mean_std[chr])):
                mean, std = bin_mean_std[chr][i]
                if std==0: continue
                zscore = (counts[chr][i]/df.loc[f, comb].sum()-mean)/std
                x += [i*50000+25000]
                y += [zscore]
            gb.scatter_on_chrom(chr, x, y, tofile=args.output+"/"+fname+"/"+chr+".html")


# Main
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""A Python implementation of NCV for NIPT. Output in tsv format.""")
    parser.add_argument('mode', default='test', type=str, help="<train/test> In train mode, -m and -f are required. In test mode, -s is required. [required, default=test]")
    parser.add_argument('-p', '--cpus', default=20, type=int, help="Number of CPUs to use for training[default=20]")
    parser.add_argument('-d', '--model', default="ncv.model", type=str, help="The model file to be write in training or read in test. [default=ncv.model]")
    parser.add_argument('-m', '--male', default="NA", type=str, help="A file that each line contains a male training sample BAM file path.")
    parser.add_argument('-f', '--female', default="NA", type=str, help="A file that each line contains a female training sample BAM file path.")
    parser.add_argument('-s', '--test', default="NA", type=str, help="A file that each line contains a test sample BAM file path.")
    parser.add_argument('-o', '--output', default=datetime.now().strftime("%Y-%m-%d_%H-%M-%S"), type=str, help="Output folder [default=<current time>]")
    args = parser.parse_args()
    if args.mode == 'test':
        if args.test == "NA":
            raise Exception("[Error] No --test detected")
        args.output = os.path.abspath(os.path.expanduser(args.output))
        args.test = os.path.abspath(os.path.expanduser(args.test))
        NCVtest(args)
    else:
        if args.male == "NA" or args.female == "NA":
            raise Exception("[Error] No --male or --female detected")
        args.male = os.path.abspath(os.path.expanduser(args.male))
        args.female = os.path.abspath(os.path.expanduser(args.female))
        NCVtrain(args)

