# run-stats.py

# program to process stats from CRAM file


import sys
import subprocess
import os
import argparse
import time
import socket
import shutil

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

