#! /usr/bin/env python
import sys, time, os, optparse, glob
import pandas as pd

def OverlapRegion(df1,df2):
	d = pd.merge(df1,df2,on="chr")
	d0 = d[(d['start_y'] >= d['start_x']) & (d['start_y'] < d['end_x']) & (d['end_y'] >= d['end_x'])].drop_duplicates()
	d0_overlap = (d0['end_x'] - d0['start_y'])/(d0['end_x']-d0['start_x'])
	d0["region"] = (d0['end_x'] - d0['start_y'])
	d0["len"] = (d0['end_x']-d0['start_x'])
	d0["coverage"] = d0_overlap
	
	d1 = d[(d['end_x'] == d['start_y'])].drop_duplicates()
	d1_overlap = 1/(d1['end_x']-d1['start_x'])
	d1["region"] = 1
	d1["len"] = (d1['end_x']-d1['start_x'])
	d1["coverage"] = d1_overlap
	
	d2 = d[(d['start_x'] >= d['start_y']) & (d['end_y'] <= d['end_x']) & (d['end_y'] > d['start_x'])].drop_duplicates()
	d2_overlap = (d2['end_y']-d2['start_x'])/(d2['end_x']-d2['start_x'])
	d2["region"] = (d2['end_y']-d2['start_x'])
	d2["len"] = (d2['end_x']-d2['start_x'])
	d2["coverage"] = d2_overlap
	
	d3 = d[(d['start_y'] > d['start_x']) & (d['end_y'] < d['end_x'])].drop_duplicates()
	d3_overlap = 1
	d3["region"] = (d3['end_y']-d3['start_y'])
	d3["len"] = (d3['end_x']-d3['start_x'])
	d3["coverage"] = d3_overlap
	
	d4 = d[(d['start_y'] < d['start_x']) & (d['end_y'] > d['end_x'])].drop_duplicates()
	d4_overlap = 1
	d4["region"] = (d4['end_x']-d4['start_x'])
	d4["len"] = (d4['end_x']-d4['start_x'])
	d4["coverage"] = d4_overlap
	
	d5 = d[(d['start_y'] < d['start_x'])].drop_duplicates()
	d5_overlap = 0
	d5["region"] = 0
	d5["len"] = (d5['end_x']-d5['start_x'])
	d5["coverage"] = d5_overlap
	
	d6 = d[(d['start_y'] > d['end_x'])].drop_duplicates()
	d6_overlap = 0
	d6["region"] = 0
	d6["len"] = (d6['end_x']-d6['start_x'])
	d6["coverage"] = d6_overlap
	
	df = d0.append(d1)
	df = d2.append(df)
	df = d3.append(df)
	df = d4.append(df)
	df = d5.append(df)
	df = d6.append(df)
	
	df0 = df.groupby(by=['chr','start_x','end_x','len']).sum().reset_index()
	del df0['start_y']
	del df0['end_y']
	df0.columns = ['chr','start','end','len','region','coverage']
	return df0


if __name__ == "__main__": 

    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)

    parser.add_option("-o","--out",dest="outfile",metavar="FILE",help="ouput overlap region file name. [default: STDOUT]",default="")
    parser.add_option("-a","--background",dest="afile",metavar="FILE",help="background file. (required)",default="")
    parser.add_option("-b","--compare",dest="bfile",metavar="FILE",help="compare file. (required)",default="")
    parser.add_option("-q","--quiet",action="store_true",dest="quiet",help="don't print progress on stderr.",default=False)

    options,args = parser.parse_args()

    if len(options.afile) == 0: parser.error("Missing background file, use -a or --background option.")
    if len(options.bfile) == 0: parser.error("Missing compare file, use -b or --compare options.")

    for i in glob.glob(options.afile):
        f1 = pd.read_table(i,sep='\t',header=None)
        df1 = pd.DataFrame(f1)
        df1.columns = ['chr','start','end']

    for j in glob.glob(options.bfile):
        f2 = pd.read_table(j,sep='\t',header=None)
        df2 = pd.DataFrame(f2)
        df2.columns = ['chr','start','end']

    out = OverlapRegion(df1,df2)

    out.to_csv(options.outfile,sep='\t',index=None,float_format='%.7f')
