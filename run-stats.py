# run-stats.py

# program to process stats from CRAM file


import sys
import subprocess
import os
import argparse
import time
import socket
import shutil
import gzip
import numpy as np

###############################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print('command failed')
        print(cmd)
        sys.exit(1)
###############################################################################
# setup paths to default programs to use and checks for required programs
def check_prog_paths(myData):        
    print('checking required programs')
    for p in ['gatk','samtools','Rscript','tabix']:
        if shutil.which(p) is None:
            s = p + ' not found in path! please fix (module load?)'
            print(s, flush=True)
            sys.exit()
        else:
            print('%s\t%s' % (p,shutil.which(p)),flush=True)
############################################################################# 
def run_flagstat(myData):
    print('\nrunning flagstat ...', flush=True)
    myData['flagStatFileName'] = myData['cramFileName'] + '.flagstat'
    
    if os.path.isfile(myData['flagStatFileName']) is False:
        cmd = 'samtools flagstat %s > %s ' % (myData['cramFileName'],myData['flagStatFileName'])
        print(cmd,flush=True)
        runCMD(cmd)
        print('Done',flush=True)
        
    
    else:
        print('%s exists, skipping flagstat' % myData['flagStatFileName'],flush=True)    
###############################################################################
def run_insertmetrics(myData):
    print('\nrunning CollectInsertSizeMetrics ...', flush=True)

    myData['sizeMetricsFileName'] = myData['cramFileName'] + '.insert_size_metrics.txt'
    myData['sizeMetricsFileHistName'] = myData['cramFileName'] + '.insert_size_metrics.hist.pdf'

    if os.path.isfile(myData['sizeMetricsFileName']) is False:
        cmd = 'gatk CollectInsertSizeMetrics -R %s -I %s -O %s -H %s' % (myData['ref'],myData['cramFileName'],
                                                                         myData['sizeMetricsFileName'],myData['sizeMetricsFileHistName'] )
        print(cmd,flush=True)
        runCMD(cmd)
        print('Done',flush=True)
    else:
        print('%s exists, skipping CollectInsertSizeMetrics' % myData['sizeMetricsFileName'],flush=True)    
###############################################################################
def run_alignmentmetrics(myData):
    print('\nrunning CollectAlignmentSummaryMetrics ...', flush=True)

    myData['alignMetricsFileName'] = myData['cramFileName'] + '.alignment_summary_metrics'

    if os.path.isfile(myData['alignMetricsFileName']) is False:
        cmd = 'gatk CollectAlignmentSummaryMetrics -R %s -I %s -O %s ' % (myData['ref'],myData['cramFileName'],
                                                                         myData['alignMetricsFileName'], )
        print(cmd,flush=True)
        runCMD(cmd)
        print('Done',flush=True)
    else:
        print('%s exists, skipping CollectAlignmentSummaryMetrics' % myData['alignMetricsFileName'],flush=True)    
###############################################################################
def run_index_cram(myData):
    print('\nrunning index cram ...', flush=True)
    craiFileName = myData['cramFileName'] + '.crai'
    if os.path.isfile(craiFileName) is False:
        cmd = 'samtools index %s' % myData['cramFileName']                                                                
        print(cmd,flush=True)
        runCMD(cmd)
        print('Done',flush=True)
    else:
        print('%s exists, skipping index cram' % craiFileName,flush=True)    
###############################################################################
def run_index_gvcf(myData):
    print('\nrunning index gvcf ...', flush=True)
    tbiFileName = myData['gvcf'] + '.tbi'

    if os.path.isfile(tbiFileName) is False:
        cmd = 'tabix -p vcf %s' %  myData['gvcf']                                                                
        print(cmd,flush=True)
        runCMD(cmd)
        print('Done',flush=True)
    else:
        print('%s exists, skipping index gvcf' % tbiFileName,flush=True)    
###############################################################################
def run_genotype_known_sites(myData):
    print('\nrunning genotype known sites ...', flush=True)
    myData['knownSitesVCF'] = myData['cramFileName'] + '.knownsites.vcf.gz'
    
    if os.path.isfile(myData['knownSitesVCF']) is False:
        cmd = 'gatk --java-options "-Xmx2g" GenotypeGVCFs '
        cmd += '-R %s -V %s -O %s ' % (myData['ref'],myData['gvcf'],myData['knownSitesVCF'] )
        cmd += ' --include-non-variant-sites --intervals %s ' % (myData['sites'])        
        print(cmd,flush=True)
        runCMD(cmd)
        print('Done',flush=True)
    else:
        print('%s exists, genotype known sites' % myData['knownSitesVCF'],flush=True)    
###############################################################################
def calc_effective_depth(myData):
    print('\nrunning calc effective depth ...', flush=True)
    myData['depthSummary'] = myData['cramFileName'] + '.knownsites.depth.txt'
    
    if os.path.isfile(myData['depthSummary']) is True:
        print('%s exists' % myData['depthSummary'],flush=True)    
        return
    autoDp = []
    xDp = []        
    inFile = gzip.open(myData['knownSitesVCF'],'rt')
    for line in inFile:
        if line[0] == '#':
            continue
        line = line.rstrip()
        line = line.split()
        
        infoField = line[7]
        infoField = infoField.split(';')
        dp = -1
        for i in infoField:
            if i[0:3] == 'DP=':
                dp = int(i.split('=')[1])
                break

        if dp == -1:
            dp = 0
        if line[0] == 'chrX':
            xDp.append(dp)
        else:
            autoDp.append(dp)
    inFile.close()
    
    outFile = open(myData['depthSummary'],'w')
        
    outFile.write('#chrom\ttotalSites\tMean\tStd\tMedian\n')
    outFile.write('Autos\t%i\t%.2f\t%.2f\t%.1f\n' % (len(autoDp),np.mean(autoDp),np.std(autoDp),np.median(autoDp) ))
    outFile.write('ChrX\t%i\t%.2f\t%.2f\t%.1f\n' % (len(xDp),np.mean(xDp),np.std(xDp),np.median(xDp) ))
    outFile.close()
###############################################################################
def summarize_stats(myData):
    myData['statsSummary'] = myData['cramFileName'] + '.stats.txt'
    
    outFile = open(myData['statsSummary'],'w')
    
    sn = myData['cramFileName'].split('/')[-1].split('.')[0]
    
    outFile.write('SampleName\t%s\n' % (sn))
    
    cramFileSize = os.path.getsize(myData['cramFileName']) / (1024*1024*1024)
    outFile.write('CramSize\t%.2f Gb\n' % cramFileSize)

    cramFileSize = os.path.getsize(myData['gvcf']) / (1024*1024*1024)
    outFile.write('GVCFSize\t%.2f Gb\n' % cramFileSize)
    
    # now get coverage stats
    
    inFile = open(myData['depthSummary'],'r')
    lines=inFile.readlines()
    inFile.close()
    
    autoLine = lines[1].rstrip().split()
    xLine = lines[2].rstrip().split()
    
    outFile.write('effectiveAutoMean\t%s\neffectiveAutoMedian\t%s\n' % (autoLine[2],autoLine[4]) )
    outFile.write('effectiveXMean\t%s\neffectiveXMedian\t%s\n' % (xLine[2],xLine[4]) )
    
    xvsAutoMean =  float(xLine[2]) / float(autoLine[2])
    xvsAutoMedian = float(xLine[4]) / float(autoLine[4])
    
    outFile.write('Mean(X/Auto)\t%.2f\n' % xvsAutoMean)
    outFile.write('Median(X/Auto)\t%.2f\n' % xvsAutoMedian)
    
    # dup mark stats
    # now get mark dup metrics
    dupMetFileName = myData['cramFileName'].replace('.cram','.sort.md.metricts.txt')
    inFile = open(dupMetFileName,'r')
    lines = inFile.readlines()
    inFile.close()

    metLine = lines[7].rstrip().split()
    dupF = metLine[8]
   
    outFile.write('fractionDup\t%s\n' % dupF) 
    
    # now get read alignment stats
    
    inFile = open(myData['alignMetricsFileName'],'r')
    lines = inFile.readlines()
    pairLine = lines[9].rstrip().split()
    inFile.close()
    
    outFile.write('totalPairedReads\t%s\n' % pairLine[1])
    outFile.write('fractionAligned\t%s\n' % pairLine[6])
    outFile.write('mismatchRate\t%s\n' % pairLine[12])
    outFile.write('indelRate\t%s\n' % pairLine[14])
    outFile.write('meanReadLen\t%s\n' % pairLine[15])
    outFile.write('fractionImproperPairs\t%s\n' % pairLine[19])
    outFile.write('fractionChimera\t%s\n' % pairLine[22])

    # now get insert len stats
    inFile = open(myData['sizeMetricsFileName'],'r')
    lines = inFile.readlines()
    inFile.close()
    
    statLine = lines[7].rstrip().split()
    
    numPairsFirstLine= int(statLine[7])
    totPairs = numPairsFirstLine
    i = 8
    while lines[i] != '\n':
       nl = lines[i].rstrip().split()
       totPairs += int(nl[7])
       i += 1
    fractionPairsAssigned = numPairsFirstLine/totPairs
       
    
    
    # check to see if


    outFile.write('pairOrientation\t%s\n' % statLine[8])
    outFile.write('fractionWithPairOrientation\t%.4f\n' % fractionPairsAssigned )    
    outFile.write('meanInsertLen\t%s\n' % statLine[5])
    outFile.write('stdInsertLen\t%s\n' % statLine[6])
    outFile.write('medianInsertLen\t%s\n' % statLine[0])
    outFile.write('madInsertLen\t%s\n' % statLine[2])
    outFile.close()  
    
    
      
    print('Summary written to',myData['statsSummary'])
###############################################################################

# SETUP

parser = argparse.ArgumentParser(description='run stats from CRAM file')

parser.add_argument('--cram', type=str,help='CRAM file',required=True)
parser.add_argument('--ref', type=str,help='genome fasta with dictionary and .fai',required=True)
parser.add_argument('--gvcf', type=str,help='sample gvcf file',required=True)
parser.add_argument('--sites', type=str,help='vcf of sites for coverage and genotyping',required=True)

args = parser.parse_args()

#####################################################################

myData = {} # dictionary for keeping and passing information
myData['cramFileName'] = args.cram
myData['ref'] = args.ref
myData['gvcf'] = args.gvcf
myData['sites'] = args.sites


# make sure programs are availble
check_prog_paths(myData)
run_flagstat(myData)
run_insertmetrics(myData)
run_alignmentmetrics(myData)
run_index_cram(myData)
run_index_gvcf(myData)
run_genotype_known_sites(myData)
calc_effective_depth(myData)
summarize_stats(myData)

