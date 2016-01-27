#!/usr/bin/env python
#Last-modified: 04 Dec 2015 12:23:40 PM

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
import time
import types
import numpy
import pandas
import matplotlib
from scipy import stats
from array import array
from ngslib import IO,DB,BedList,Utils
from matplotlib.pyplot import plot,subplots,savefig

# ------------------------------------
# constants
# ------------------------------------

Chroms = {'hg19': "chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrM chrX chrY".split(),
         }

gversion = 'hg19'

# Seed related variables
bitsToBases = {0:'T', 1:'C', 2:'A', 3:'G'}
byteTable = {}
twoBytesTable = {}
seedsTable = {}
bases = 'ACTG'

# dumpfile header
header_type = numpy.uint64
header_length = 24

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

class Algorithms(object):
    def aRRA(rankedlst,geneset):
        '''
        '''
        M = float(len(rankedlst))
        genes = {gid:i/M for i,gid in enumerate(rankedlst)}
    aRRA=staticmethod(aRRA)

class SUtils(object):
    def fetchSeq(genome,gannos,utrfile,lncfile):
        '''
        Fetch sequence from genome according to transcriptome annotation.
        Parameters:
            genome: string
                genome in fasta format
            gannos: string or list
                transcriptiome annotation(s) in genepred format
            utrfile: string
                output file for 3'UTR
            lncfile: string
                output file for lncRNAs
        '''
        gdb = DB(genome,'fasta')
        utrfh = open(utrfile,'w')
        lncfh = open(lncfile,'w')
        chroms = Chroms[gversion]
        if isinstance(gannos,str):
            gannos = [gannos]
        for ganno in gannos:  
            for gene in IO.BioReader(ganno,'genepred'):
                if gene.chrom not in chroms:
                    continue
                gid = gene.id if gene.proteinid == "" else gene.proteinid
                ofh = lncfh # by default it is lncRNA. NOTE: we didn't check seq length.
                if gene.txstart != gene.txstop: # protein coding 
                    gene = gene.getUTR3()
                    if not gene:
                        continue
                    ofh = utrfh
                seq = gene.fetchDB(gdb)
                print >>ofh, ">{0}:{1}\n{2}".format(gid,gene.id,seq)
        gdb.close()
        utrfh.close()
        lncfh.close()
    fetchSeq=staticmethod(fetchSeq)
    def getSeeds(fafile,sense=True,start_at=2,stop_at=8):
        if sense:
            for fa in IO.BioReader(fafile,'fasta'):
                yield fa.id,fa.seq.seq[start_at-1:stop_at-1]
        else:
            for fa in IO.BioReader(fafile,'fasta'):
                yield fa.id,Utils.rc(fa.seq.seq)[start_at-1:stop_at-1]
    def touchtime(s='# {0}:'):
        return s.format(time.ctime())
    touchtime=staticmethod(touchtime)

class SeedUtils(object):
    '''
    Seeds related functions.
    '''
    scnts = None # seed counts given a siRNA file.
    def createTables():
        ''' Create byteTable and twoBytesTables. '''
        for x in xrange(2**8):
            c = (x >> 4) & 0xf
            f = x & 0xf
            cc = (c >> 2) & 0x3
            cf = c & 0x3
            fc = (f >> 2) & 0x3
            ff = f & 0x3
            s = ''.join(map(lambda x:bitsToBases[x], (cc, cf, fc, ff)))
            byteTable[x] = s
        for x in xrange(2**16):
            c = (x >> 8) & 0xff
            f = x & 0xff
            if x >>14 !=0: # if not 'T', set it to 1
                x = (x & 0x3fff) + 0x4000
            twoBytesTable[byteTable[c] + byteTable[f]] = x
            if x< 16384: # 2**14
                seedsTable[x] = byteTable[c][1:] + byteTable[f]
    createTables=staticmethod(createTables)    
    def parseHeader(fh):
        '''
        '''
        header = numpy.frombuffer(fh.read(header_length),dtype=header_type,count=3)
        gids = fh.read(int(header[2]))
        genes = []
        lengths = []
        for item in gids.split(';'):
            idx = item.rfind(':')
            genes.append(item[:idx])
            lengths.append(int(item[idx+1:]))
        return genes,lengths
    parseHeader=staticmethod(parseHeader)
    def countSeedInSiRNA(siRNASeeds):
        '''
        Read si(h)RNA file in fasta format. Calculate the frequency of each seed.
        '''
        scnts = numpy.zeros(2**14,dtype=numpy.uint16)
        for seed in siRNASeeds:
            scnts[twoBytesTable['T'+seed]] += 1
        return scnts
    countSeedInSiRNA=staticmethod(countSeedInSiRNA)
    def findSeeds(seq,sary):    
        sary.fill(0)
        seq = Utils.rc(seq)
        for i in range(len(seq)-7):
            sary[twoBytesTable[seq[i:i+8]]] += 1
    findSeeds=staticmethod(findSeeds)
    def parseSeeds(mat,seeds):
        '''
        Calculate 4 types of seed matches.
        '''
        N = mat.shape[0] # number of genes
        if N: # N!= 0        
            for x in range(2**14):
                seed = seedsTable[x]
                seed, P8 = seed[:6], seed[-1]
                NP8 = [c for c in bases if c != P8]
                seeds[:,0,x] = mat[:,twoBytesTable['T'+seed+P8]]
                seeds[:,1,x] += mat[:,twoBytesTable['A'+seed+P8]] # same index for 'ACG' at first position
                for idx in [twoBytesTable['T'+seed+c] for c in NP8]:
                    seeds[:,2,x] += mat[:,idx]
                for c in NP8:
                    seeds[:,3,x] += mat[:,twoBytesTable['A'+seed+c]] # same index for 'ACG' at first position
    parseSeeds=staticmethod(parseSeeds)
    def parseGeneSeeds(gfile,dumpfile,N=10):
        '''
        Parse seeds on gene sequences.
        Parameters:
            gfile: string
                Gene sequence file in Fasta format. 
            N: int
                Parse N sequence at a time instead all to reduce memory usage. [Default=10000]
                N = 10000, memory usage for the two huge matrix is up to 1.83G.
        '''
        # create seed tables
        SeedUtils.createTables() 
        # read genes
        ofh = gzip.open(dumpfile,'wb')
        genes = [fa for fa in IO.BioReader(gfile,'fasta')]
        # write headers
        gids = ';'.join(["{0}:{1}".format(fa.id,len(fa)) for fa in genes])
        header = numpy.array([0x19840405,len(genes),len(gids)],dtype=numpy.uint64)
        ofh.write(header.data)
        ofh.write(gids)
        # parse genes 
        cnt = 0
        mat = numpy.zeros((N,2**15),dtype=numpy.uint16)
        seeds = numpy.zeros((N,4,2**14),dtype=numpy.uint16) # N x type x seeds
        for fa in genes:
            SeedUtils.findSeeds(fa.seq.seq.upper(),mat[cnt%N])
            cnt += 1
            if cnt%(N)==0:
                # calculate seeds for the N genes
                SeedUtils.parseSeeds(mat,seeds)
                ofh.write(seeds.data)
                print >>sys.stderr, touchtime(), "Parsed {0} genes .".format(cnt)
                seeds.fill(0)
                mat.fill(0)
        # calculate seeds for the rest of genes
        rest = cnt%N
        if rest:
            SeedUtils.parseSeeds(mat[:rest],seeds[:rest])
            ofh.write(seeds[:rest].data)
        print >>sys.stderr, Utils.touchtime(), "Parsed {0} genes.       ".format(cnt)
        ofh.close()
    parseGeneSeeds=staticmethod(parseGeneSeeds)
    def calOffTargetEffect(gdumpfile,ofedumpfile,weight=[1.,1.,1.,1.],N=10):
        '''
        Calculate the off-target effect of between seeds and genes.
        In the final version, this will be merged in the parseGeneSeeds function.
        '''
        with gzip.open(gdumpfile,'rb') as fh:
            ofh = gzip.open(ofedumpfile,'wb')
            header = numpy.frombuffer(fh.read(header_length),dtype=header_type,count=3)
            num_gene = int(header[1]) # number of genes
            # check signature
            if header[0] == 0x19840405:
                gids = fh.read(int(header[2]))
                cnt = 0
                ofe_mat = numpy.zeros((N,2**14),dtype=numpy.float16) # N*seeds
                while cnt<num_gene:
                    if cnt+N>num_gene:
                        N = num_gene%N
                        if N == 0:
                            break
                    # calculate OffTargetEffect
                    seeds = numpy.frombuffer(fh.read(2**17*N),dtype=numpy.uint16,count=2**16*N) # 2**14*4*2*N bytes
                    seeds.reshape((N,4,2**14))
                    for i in range(4):
                        ofe_mat[:N,:] += weight[i]*seeds[:N,i,:]
                        ofh.write(ofe.data)
            ofh.close()
    calOffTargetEffect=staticmethod(calOffTargetEffect)
    def aRRA(siRNASeeds,geneSeeds,a = 0.05):
        '''
        aRRA run in parallel on matrix.
        Parameters:
            siRNASeeds: list of strings
                siRNA seed sequences in order of phenotype.
            geneSeeds: string
                Dumped gene seeds file.
            a: float
                Consider only top a percent of siRNAs. 
        '''
        # Count seeds in siRNA file
        if SeedUtils.scnts is None:
            SeedUtils.scnts = SeedUtils.countSeedInSiRNA(siRNASeeds)
        # Count total number of siRNAs linked to each gene
        tcnts = numpy.sum(numpy.multiply(geneSeeds!=0,SeedUtils.scnts),axis=1)
        # number of siRNAs
        N = sum(SeedUtils.scnts)
        # current seed counts
        ccnt = numpy.zeros(2**14,dtype=numpy.uint16)
        # pvalues preset to 1.0
        pvalues = numpy.ones(geneSeeds.shape[0],dtype=float)
        # current siRNA-gene links
        ck = numpy.zeros(geneSeeds.shape[0],dtype=numpy.uint16)
        # top alpha of the list
        aN = int(a*N)
        for i,seed in zip(range(aN),siRNASeeds):
            x = twoBytesTable['T'+seed ]   # seed index
            ccnt[x] += 1                                # add seed to ccnt
            idx = geneSeeds[:,x] !=0                    # indices of genes having this seed
            ck[idx] += 1                                # add the siRNA-gene link to ck
            pvalues[idx] = numpy.min([pvalues[idx],stats.beta.cdf((i+1.)/N,ck[idx],tcnts[idx]+1-ck[idx])],axis=0)
        return -numpy.log(pvalues)
    aRRA=staticmethod(aRRA)
    def evdplot(df,outfile):
        '''
        Distribution of escores.
        '''
        matplotlib.use('pdf')
        f, (ax1, ax2, ax3) = subplots(1,3,figsize=[10,5])
        # s score distribution
        df.escore.hist(ax=ax1,normed=1,bins=50,histtype='stepfilled', alpha=0.2,label='EScore')
        c,loc,scale = stats.genextreme.fit(df.escore)
        rv = stats.genextreme(c,loc=loc,scale=scale)
        x = numpy.linspace(0,df.escore.max(),50)
        ax1.plot(x, rv.pdf(x), 'k-', lw=1, label='EVD')
        ax1.legend(loc='best')
        ax1.set_xlabel('Enrichment score')
        ax1.set_ylabel('Probability')
    
        # cummulative distribution
        df.escore.hist(ax=ax2,cumulative=True, normed=1, bins=50,histtype='stepfilled', alpha=0.2,label='EScore')
        ax2.plot(x, rv.cdf(x), 'k-', lw=1, label='EVD')
        ax2.legend(loc='best')
        ax2.set_xlabel('Enrichment score')
        ax2.set_ylabel('Probability')
    
        # fitness test
        nbin = int(round(1+numpy.log2(df.escore.size)))
        x = numpy.linspace(0,df.escore.max(),nbin+1)
        y = rv.cdf(x)
        counts, bin_edges = numpy.histogram(df.escore, bins=nbin, normed=False)
        counts = [df.escore.size-len(df.escore.nonzero()[0])] + list(counts)
        cdf = numpy.cumsum(counts)
        cdf = cdf/float(max(cdf))
    
        kst,ksp = stats.ks_2samp(y,cdf)
        chit,chip = stats.chisquare(cdf,y)
        ax3.plot(bin_edges,y,'r-',label='EVD')
        ax3.plot(bin_edges, cdf,'b-',label='EScore')
        ax3.legend(loc='best')
        ax3.text(df.escore.max()/4,0.3,"KS test:\nstat={0:.2f},p={1:.2e}\nChiSquare test:\nstat={2:.2f},p={3:.2e}".format(kst,ksp,chit,chip))
        ax3.set_xlabel('Enrichment score')
        ax3.set_ylabel('Probability')

        savefig(outfile,format='pdf')
        return rv
    evdplot=staticmethod(evdplot)
    def calPvalues(siRNASeeds,seedfile,outfile,method='RRA',N=10000):
    
        # read headers
        print >>sys.stderr, Utils.touchtime(),"reading header"
        fh = gzip.open(seedfile,'rb')
        genes,lengths = SeedUtils.parseHeader(fh)
        num_genes = len(genes)
    
        # aRRA
        print >>sys.stderr, Utils.touchtime(),"start analyzing data."
        method = {'RRA':SeedUtils.aRRA}.get(method,SeedUtils.aRRA)
        escores = numpy.zeros(num_genes)
        for i in range(0,num_genes,N):
            n = N if i+N < num_genes else num_genes%N
            if n:
                geneSeeds = numpy.sum(numpy.frombuffer(fh.read(2**17*n),dtype=numpy.uint16,count=2**16*n).reshape((n,4,2**14)),axis=1)
                escores[i:i+n] = method(siRNASeeds,geneSeeds)
                print Utils.touchtime(), "analyzed", i+n
        fh.close()
    
        # Writting result to file
        print >>sys.stderr, Utils.touchtime(), 'writing result to file.'
        df = pandas.DataFrame({'gid':genes,'length':lengths, 'escore':escores})
        rv = SeedUtils.evdplot(df,outfile+".pdf")
        df.loc[:,'pvalue'] =  rv.pdf(df.escore)
        df = df.sort_values(by='pvalue',ascending=False)
        df.to_csv(outfile,index=False,sep='\t',columns=['gid','length','escore','pvalue'])
        print >>sys.stderr, Utils.touchtime(), "finished."
        
        return df
    calPvalues=staticmethod(calPvalues)
    
# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        sys.exit("Example:"+sys.argv[0]+" file1 file2... ")
    for item in IO.BioReader(sys.argv[1],ftype='bed'):
        print item

