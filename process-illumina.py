# process-illumina.py

# program to run the alignment of Illumina WGS data
# goes from fastq.gz to CRAM + GVCF 

import dogmap
import os
import sys
import argparse
import time

# SETUP

parser = argparse.ArgumentParser(description='process Illumina Reads')

parser.add_argument('--sn', type=str,help='sample name',required=True)
parser.add_argument('--lib', type=str,help='library name',required=True)

parser.add_argument('--fq1', type=str,help='fastq of read1',required=True)
parser.add_argument('--fq2', type=str,help='fastq of read1',required=True)

parser.add_argument('--ref', type=str,help='genome fasta with dictionary and .fai',required=True)
parser.add_argument('--refBWA', type=str,help='genome fasta BWA-MEM2 index',required=True)

parser.add_argument('-t', type=int,help='number of threads to use',required=True)

parser.add_argument('--tmpdir', type=str,help='tmp dir for running',required=True)
parser.add_argument('--finaldir', type=str,help='final dir for output',required=True)

parser.add_argument('--knownsites', type=str,help='vcf of known sites for BQSR',required=True)


args = parser.parse_args()

#####################################################################

myData = {} # dictionary for keeping and passing information
myData['sampleName'] = args.sn
myData['libName'] = args.lib
myData['fq1'] = args.fq1
myData['fq2'] = args.fq2
myData['tmpDir'] = args.tmpdir
myData['finalDir'] = args.finaldir
myData['ref'] = args.ref
myData['refBWA'] = args.refBWA
myData['threads'] = args.t
myData['knownSitesVCF'] = args.knownsites


if myData['tmpDir'][-1] != '/':
   myData['tmpDir'] += '/'
if myData['finalDir'][-1] != '/':
   myData['finalDir'] += '/'


myData['logFileName'] = myData['finalDir'] + myData['sampleName'] + '.map.log'
myData['completeToken'] = myData['finalDir'] + myData['sampleName'] + '.map.complete'


# if log file exists, then there is partial processing so not sure we want to redo and overwrite
# safe to just quite and letter user deal with it

if os.path.isfile(myData['logFileName']) is True:
    print('ERROR!!!')
    print('%s exists.  Do you really want to rerun this pipeline?' % myData['logFileName'] )
    print('ERROR!!!')
    sys.exit()

myData['logFile'] = open(myData['logFileName'],'w')

# add initial infoto log
dogmap.init_log(myData)

# make sure programs are availble
dogmap.check_prog_paths(myData)
dogmap.check_dir_space(myData)
dogmap.run_bwa_mem2(myData)
dogmap.run_mdspark(myData)

dogmap.run_bqsr(myData)
dogmap.run_haplotypecaller(myData)

# clean up and get elapsed time!
dogmap.remove_tmp_dir(myData)

cmd = 'touch %s ' % myData['completeToken']
print(cmd,flush=True)
myData['logFile'].write(cmd + '\n')
dogmap.runCMD(cmd)

myData['endTime'] = time.localtime()
myData['tEnd'] = time.time()
t = time.strftime("%a, %d %b %Y %H:%M:%S", myData['endTime'])
myData['logFile'].write('\nEnded!\n%s\n' % t)

elapsedTime = myData['tEnd'] - myData['tStart'] # this is in nanoseconds??

# convert to minutes
elapsedTime= elapsedTime / 60
# convert to hours
elapsedTime = elapsedTime / 60

#t = time.strftime("%H:%M:%S",elapsedTime )
myData['logFile'].write('Elapsed time:\n%s hours\n' % elapsedTime)
myData['logFile'].close()



