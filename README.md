# A NCV implementation by Luyu

### Installation
1. Install **BWA** and **samtools**
2. Prepare genome file **hg38.fa**
    1. Download `wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`
    2. Unzip `gunzip hg38.fa.gz`
    3. Build index `bwa index hg38.fa` 
3. Install `python2`/`python3`
5. Install `numpy`, `pandas`, `pysam`, `tqdm`, `pp` by pip
6. Download all files here

## Structure

`NCV_ver_1.6.py`: main executable python script.
`genome_browser.py`: auxiliary file to generate chromosome view html.
`batch3_females` and `batch3_males`: example of the format for `-f`/`-m`/`-s` arguments.
`ncv.model`: a pre-trained model to quick start.
`hg38_cytoBand.txt`: the hg38 reference cytoband annotation file.

### Usage
1. use BWA to map raw reads to reference genome  and sort them by samtools `bwa mem -t 20 hg38.fa sample.fq.gz | samtools sort -@ 20 - -o sample.bam`
2. index BAM fiels `samtools index sample.bam` 
3. train or test
    1. train the model with Normal NIPT samples:
`python NCV_1.5.py train -m male.lst -f female.lst`
    2. test a batch of NIPT samples:
`python NCV_1.5.py test -s test.lst`
 
---
Type `python NCV_ver_1.6.py --help` to show help information:
```
usage: NCV_ver_1.6.py [-h] [-p CPUS] [-d MODEL] [-m MALE] [-f FEMALE]
                      [-s TEST] [-o OUTPUT]
                      mode

A Python implementation of NCV for NIPT. Output in tsv format.

positional arguments:
  mode                  <train/test> In train mode, -m and -f are required. In
                        test mode, -s is required. [required, default=test]

optional arguments:
  -h, --help            show this help message and exit
  -p CPUS, --cpus CPUS  Number of CPUs to use for training[default=20]
  -d MODEL, --model MODEL
                        The model file to be write in training or read in
                        test. [default=ncv.model]
  -m MALE, --male MALE  A file that each line contains a male training sample
                        path.
  -f FEMALE, --female FEMALE
                        A file that each line contains a female training
                        sample path.
  -s TEST, --test TEST  A file that each line contains a test sample path.
  -o OUTPUT, --output OUTPUT
                        Output filename prefix [default=NCV]
```
