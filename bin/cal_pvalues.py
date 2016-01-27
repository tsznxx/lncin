#!/usr/bin/env python
#Last-modified: 02 Dec 2015 03:36:17 PM

#         Module/Scripts Description
# 
# Copyright (c) 2008 Yunfei Wang <Yunfei.Wang1@utdallas.edu>
# 
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @status:  experimental
# @version: 1.0.0
# @author:  Yunfei Wang
# @contact: yfwang0405@gmail.com

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import gzip
import numpy
import pandas
import argparse
from scipy import stats
from lncin import SeedUtils,Utils
from ngslib import IO
# ------------------------------------
# constants
# ------------------------------------

TRCfile = "/scratch/bcb/ywang52/TData/siRNA-offtarget/TRC_public_05April.fa"
gfile = "/scratch/bcb/ywang52/TData/siRNA-offtarget/020prediction/miranda/gencodev19_mRNA_UTR3.fa"
dumpfile = "/scratch/bcb/ywang52/TData/siRNA-offtarget/020prediction/TargetScan/TargetScanS/gencodev19_mRNA_UTR3.dump.gz"


# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Calculate enrichment scores and p values of given an ordered siRNA list.',epilog='dependency lncin')
    p._optionals.title = "Options"
    p.add_argument("-i",dest="siRNA",type=str,metavar="siRNA.fa", required=True, help="Ordered siRNA file in fasta format.")
    p.add_argument("-m",dest="method",type=str,metavar="method", choices=['RRA'],default='RRA',help="Method used to calculate the enrichment score. [default='RRA']")
    p.add_argument("-s",dest="seedfile",type=str,metavar='seedfile',default=dumpfile,help="Dumped seed file. [default='gencodev19_mRNA_UTR3.dump.gz']")
    p.add_argument("-a",dest="alpha", type=float,metavar='0.05' ,default=0.05,help="Top percentage of siRNAs considered in calculation. [Default=0.05].")
    p.add_argument("-o",dest="outfile",type=str,metavar="outfile",default="stdout",help="Output file. [default= 'stdout'].")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

# ------------------------------------
# Classes
# ------------------------------------


# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    args = argParser()
    
    # Initiation
    SeedUtils.createTables()
    siRNASeeds = [line.split()[1][1:8] for line in open(args.siRNA)]

    df = SeedUtils.calPvalues(siRNASeeds,dumpfile,args.outfile)

