import plotly.graph_objs as go
from plotly.subplots import make_subplots
import pandas as pd
import plotly.colors

class GenomeBroswer:
    def __init__(self, cyto_filepath="hg38_cytoBand.txt"):
        self.cytoband_df = pd.read_table(cyto_filepath, names=["chrom", "start", "end", "label", "cyto"])
        cyto_abbr = ['gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100', 'acen', 'gvar', 'stalk']
        self.p_cyto = {cyto_abbr[i]: plotly.colors.sequential.Blues[i+1] for i in range(8)}
        self.q_cyto = {cyto_abbr[i]: plotly.colors.sequential.Purples[i+1] for i in range(8)}

    def fill_banner_ideogram(self, chrom, fig):
        chrom_cyto_df = self.cytoband_df[self.cytoband_df['chrom']==chrom]
        for _, record in chrom_cyto_df.iterrows():
            fig.add_trace(go.Scatter(x=[record.start, record.end, record.end, record.start], 
                                            y=[1,1,-1,-1], 
                                            fill="toself", 
                                            fillcolor=self.p_cyto[record.cyto] if record.label[0]=='p' else self.q_cyto[record.cyto],
                                            showlegend=False,
                                            name=record.label+"("+record.cyto+")",
                                            mode="none"), row=1, col=1)

    def scatter_on_chrom(self, chrom, x, y, tofile=None):
        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.1, 0.9], vertical_spacing=0)
        self.fill_banner_ideogram(chrom, fig)
        color = ["red" if abs(i)>3 else "black" for i in y]
        fig.add_trace(go.Scatter(x=x, y=y, marker={"color": color}, mode="markers", name=chrom, showlegend=False), row=2, col=1)
        fig.update_layout(xaxis={"title":chrom,
                                 "side":"top",
                                 "showticklabels":True,
                                 "range":[0,self.cytoband_df[self.cytoband_df['chrom']==chrom]['end'].max()]},
                          xaxis2={"showticklabels": False},
                          yaxis={"range":[-1,1],
                                "showticklabels":False})
        if tofile: 
            fig.write_html(tofile)
        else: 
            return fig
