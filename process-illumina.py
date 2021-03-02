# process-illumina.py

# program to run the alignment of Illumina WGS data
# goes from fastq.gz to CRAM + GVCF 

import sys
import subprocess
import os
import argparse
import time
import socket
import shutil

###############################################################################
def read_fasta_file_to_dict(fastaFile):
    myDict = {}
    inFile = open(fastaFile,'r')
    line = inFile.readline()
    line = line.rstrip()
    if line[0] != '>':
        print('ERROR, FILE DOESNNOT START WITH >')
        sys.exit()
    myName = line[1:].split()[0]  # ignore other info
    myDict[myName] = {}
    myDict[myName]['seq'] = ''
    myDict[myName]['seqLen'] = 0    
    mySeq = ''
    while True:
        line = inFile.readline()
        if line == '':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            break
        line = line.rstrip()
        if line[0] == '>':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            myName = line[1:]
            myDict[myName] = {}
            myDict[myName]['seq'] = ''
            myDict[myName]['seqLen'] = 0    
            mySeq = ''
            continue
        mySeq += line
    inFile.close()
    return myDict
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
# Helper function to run commands, handle return values and print to log file
def runCMD_output(cmd):
    val = subprocess.Popen(cmd, text=True, shell=True, stdout = subprocess.PIPE)
    resLines = []
    for i in val.stdout:
       i = i.rstrip()
       resLines.append(i)
    return resLines
#############################################################################        
#############################################################################        
# Helper function to read in information from genome .fai file and return
# a dictionary containing chrom names and lengths
def read_chrom_len(faiFileName):
    chromLens = {}
    inFile = open(faiFileName,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        chromLens[line[0]] = int(line[1])
    inFile.close()
    return chromLens    
############################################################################# 
# setup paths to default programs to use and checks for required programs
def check_prog_paths(myData):    
    
    myData['logFile'].write('\nChecking for required programs...\n')
    
    for p in ['bwa-mem2','gatk','samtools','parallel']:
        if shutil.which(p) is None:
            s = p + ' not found in path! please fix (module load?)'
            print(s, flush=True)
            myData['logFile'].write(s + '\n')
            myData['logFile'].close()        
            sys.exit()
        else:
            myData['logFile'].write('%s\t%s\n' % (p,shutil.which(p)))
            
    myData['logFile'].flush()              
############################################################################# 
def init_log(myData):
    k = list(myData.keys())
    k.sort()
    

    myData['startTime'] = time.localtime()
    myData['tStart'] = time.time()
    t = time.strftime("%a, %d %b %Y %H:%M:%S", myData['startTime'])
        
    myData['logFile'].write(t + '\n')
    
    hn = socket.gethostname()
    myData['logFile'].write('Host name: %s\n' % hn)
    myData['logFile'].write('\nInput options:\n')
    for i in k:
        if i in ['logFile']:
            continue        
        myData['logFile'].write('%s\t%s\n' % (i,myData[i]))
        
        
    myData['logFile'].flush()  
############################################################################# 
def check_dir_space(myData):
    myData['logFile'].write('\nchecking file systems\n')
    
    # check tmp dir
    if os.path.isdir(myData['tmpDir']) is False:
        s = myData['tmpDir'] + ' is not found! please check'
        print(s,flush=True)
        myData['logFile'].write(s + '\n')
        myData['logFile'].close()        
        sys.exit()

    if os.path.isdir(myData['finalDir']) is False:
        s = myData['finalDir'] + ' is not found! please check'
        print(s,flush=True)
        myData['logFile'].write(s + '\n')
        myData['logFile'].close()        
        sys.exit()
        
    cmd = 'df -h %s' % myData['tmpDir']
    o = runCMD_output(cmd)
    myData['logFile'].write(cmd + '\n')
    myData['logFile'].write(o[0] + '\n')
    myData['logFile'].write(o[1] + '\n')
    myData['logFile'].write('\n')

    cmd = 'df -h %s' % myData['finalDir']
    o = runCMD_output(cmd)
    myData['logFile'].write(cmd + '\n')
    myData['logFile'].write(o[0] + '\n')
    myData['logFile'].write(o[1] + '\n')
    myData['logFile'].flush()  
    

############################################################################# 
def run_bwa_mem2(myData,run=True):
# run bwa mem2, only option is to not run program, use in reset mode
    s = 'Starting bwa mem'
    print(s,flush=True)
    myData['logFile'].write('\n' + s + '\n')
    myData['logFile'].flush()              
    
    myData['workingBaseDir'] = myData['tmpDir'] + myData['sampleName']
    if os.path.isdir(myData['workingBaseDir']) is True:
        s = '%s exists!' % (myData['workingBaseDir'])
        print(s,flush=True)
        myData['logFile'].write(s + '\n')
    else:
       cmd = 'mkdir %s' % myData['workingBaseDir']
       myData['logFile'].write(cmd + '\n')
       runCMD(cmd)
    
    myData['workingBaseDir'] += '/'   
    myData['bwaMEMBam'] = myData['workingBaseDir'] + myData['sampleName'] + '.bam'
    
    cmd = 'bwa-mem2 mem -K 100000000  -t %i ' % (myData['threads'])
    rg = '\'@RG\\tID:%s\\tSM:%s\\tLB:%s\\tPL:ILLUMINA\'' % (myData['sampleName']+'_'+myData['libName'],myData['sampleName'],myData['libName'])
    cmd += ' -R %s ' % rg
    cmd += '%s %s %s' % (myData['refBWA'], myData['fq1'],myData['fq2'])
    cmd += ' | samtools view -bS - >  %s ' % myData['bwaMEMBam']
    
    if run is True:    
        print(cmd)
        myData['logFile'].write(cmd + '\n')
        runCMD(cmd)        
    else:
        s = 'skipping run_bwa_mem2'
        print(s,flush=True)
        myData['logFile'].write(s + '\n')    
        
    myData['logFile'].flush()              
############################################################################# 
def run_mdspark(myData,run=True):
    s = 'Starting run_mdspark'
    print(s,flush=True)
    myData['logFile'].write('\n' + s + '\n')
    myData['logFile'].flush()              
    
    myData['MDbam'] = myData['workingBaseDir'] + myData['sampleName'] + '.sort.md.bam'
    myData['dupMetsFile'] = myData['finalDir'] + myData['sampleName'] + '.sort.md.metricts.txt'
    
    
    cmd = 'gatk MarkDuplicatesSpark -I %s -O %s -M %s ' % (myData['bwaMEMBam'],myData['MDbam'],myData['dupMetsFile'])
    cmd += ' --tmp-dir %s ' % myData['tmpDir']
    cmd += ' --conf ''spark.executor.cores=%i'' ' % myData['threads']
    cmd += ' --conf ''spark.local.dir=%s'' ' % myData['tmpDir']

    if run is True:    
        print(cmd)
        myData['logFile'].write(cmd + '\n')
        myData['logFile'].flush()        
        runCMD(cmd)        
    else:
        s = 'skipping run_bwa_mem2'
        print(s,flush=True)
        myData['logFile'].write(s + '\n')    
        
    myData['logFile'].flush()              
############################################################################# 
def run_bqsr(myData,run=True):
    #setup, run, and apply BQSR
    s = 'starting run_bqsr'
    print(s,flush=True)
    myData['logFile'].write('\n' + s + '\n')
    myData['logFile'].flush()              
   


    # first, need to make list of intervals to procsess!
    setup_intervals(myData)
    
    # now, can write out cmd for each intervals
    
    # make intervals files
    myData['bqsrDataDir'] = myData['workingBaseDir'] + 'bqsr-recal'
    
    if os.path.isdir(myData['bqsrDataDir']) is True:
        s = '%s exists!' % (myData['bqsrDataDir'])
        print(s,flush=True)
        myData['logFile'].write(s + '\n')
    else:
       cmd = 'mkdir %s' % myData['bqsrDataDir']
       myData['logFile'].write(cmd + '\n')
       runCMD(cmd)    
    myData['bqsrDataDir'] += '/'
    
    myData['bqsrJobsFileName'] = myData['workingBaseDir'] + 'baserecal.jobs.txt'
    outFile = open(myData['bqsrJobsFileName'],'w')
    
    listOfBQSRRecalFiles = []
    
    for intervalFileName in myData['bqsrIntervalsFiles']:
        cName = intervalFileName.split('/')[-1].split('.')[0]
        outDataName = myData['bqsrDataDir'] + cName + '.recal_data.table'
        
        cmd = 'gatk --java-options "-Xmx4G" BaseRecalibrator '
        cmd += ' --tmp-dir %s ' % myData['tmpDir']
        cmd += ' -I %s ' % myData['MDbam']
        cmd += ' -R %s' % myData['ref']
        cmd += ' --intervals %s ' % intervalFileName
        cmd += ' --known-sites %s ' % myData['knownSitesVCF']
        cmd += ' -O %s ' % outDataName        
        cmd += '\n'
        outFile.write(cmd)
        
        listOfBQSRRecalFiles.append(outDataName)
        
    outFile.close()
    
    s = 'list of BaseRecalibrator written to %s' % myData['bqsrJobsFileName']
    print(s,flush=True)
    myData['logFile'].write('\n' + s + '\n')
    myData['logFile'].flush()
    
    
    cmd = 'parallel --jobs %i  < %s' % (myData['threads'],myData['bqsrJobsFileName'])              

    if run is True:    
        print(cmd)
        myData['logFile'].write(cmd + '\n')        
        runCMD(cmd)        
        myData['logFile'].flush()
    else:
        s = 'skipping parallel BaseRecalibrator'
        print(s,flush=True)
        myData['logFile'].write(s + '\n')    
        myData['logFile'].flush()

    
    # next, need to gather up the BaseRecalibrator reports
    myData['bqsrReportsFileList'] = myData['bqsrJobsFileName'] + '.reports.list'
    myData['bqsrReportsGatheredFile'] = myData['workingBaseDir']  + 'bqsrgathered.reports.list'
    outFile = open( myData['bqsrReportsFileList'],'w')
    for f in listOfBQSRRecalFiles:
        outFile.write('%s\n' % f)    
    outFile.close()
    
    cmd = 'gatk --java-options "-Xmx6G" GatherBQSRReports '
    cmd += ' --input %s ' % myData['bqsrReportsFileList']
    cmd += ' --output %s ' % myData['bqsrReportsGatheredFile']
    
    if run is True:    
        print(cmd)
        myData['logFile'].write(cmd + '\n')        
        runCMD(cmd)        
        myData['logFile'].flush()
    else:
        s = 'skipping GatherBQSRReports'
        print(s,flush=True)
        myData['logFile'].write(s + '\n')    
        myData['logFile'].flush()
    
    # now, apply BQSR by chromosome again...
    myData['ApplyBQSRJobsFileName'] = myData['workingBaseDir'] + 'applyBQSR.jobs.txt'
    outFile = open(myData['ApplyBQSRJobsFileName'],'w')
    
    for intervalFileName in myData['bqsrIntervalsFiles']:
        cName = intervalFileName.split('/')[-1].split('.')[0]
        outBAMName = myData['bqsrDataDir'] + cName + '.bqsr.cram'        
        cmd = 'gatk --java-options "-Xmx4G" ApplyBQSR '
        cmd += ' --tmp-dir %s ' % myData['tmpDir']
        cmd += ' -I %s ' % myData['MDbam']
        cmd += ' -R %s' % myData['ref']
        cmd += ' -O %s ' % outBAMName        
        cmd += ' --intervals %s ' % intervalFileName
        cmd += ' --bqsr-recal-file %s ' % myData['bqsrReportsGatheredFile']
        cmd += ' --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20  --static-quantized-quals 30  '
        
        cmd += '\n'
        outFile.write(cmd)
    
    # do unmapped
    outBAMName = myData['bqsrDataDir'] + 'unmapped' + '.bqsr.cram'        
    cmd = 'gatk --java-options "-Xmx4G" ApplyBQSR '
    cmd += ' --tmp-dir %s ' % myData['tmpDir']
    cmd += ' -I %s ' % myData['MDbam']
    cmd += ' -R %s' % myData['ref']
    cmd += ' -O %s ' % outBAMName        
    cmd += ' --intervals %s ' % 'unmapped'
    cmd += ' --bqsr-recal-file %s ' % myData['bqsrReportsGatheredFile']
    cmd += ' --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20  --static-quantized-quals 30  '

    cmd += '\n'
    outFile.write(cmd)
    
    
    outFile.close()
    
    s = 'list of ApplyBQSRJobs written to %s' % myData['ApplyBQSRJobsFileName']
    print(s,flush=True)
    myData['logFile'].write('\n' + s + '\n')
    myData['logFile'].flush()
    
    cmd = 'parallel --jobs %i  < %s' % (myData['threads'],myData['ApplyBQSRJobsFileName'])              

    if True is True:    
        print(cmd)
        myData['logFile'].write(cmd + '\n')        
        runCMD(cmd)        
        myData['logFile'].flush()
    else:
        s = 'skipping parallel ApplyBQSRJobsFileName'
        print(s,flush=True)
        myData['logFile'].write(s + '\n')    
        myData['logFile'].flush()
     
    # next, need to merge together CRAM files

    # GatherBamFiles, need them in order to be used
    
    
    
    
    
    
############################################################################# 
############################################################################# 
def setup_intervals(myData):
    # read in contig names
    contigNames = []
    inFile = open(myData['ref'] + '.fai')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        c = line[0]
        contigNames.append(c)
    inFile.close()
    s = 'read in %i chrom names, working on assiging to intervals' % len(contigNames)
    print(s,flush=True)
    myData['logFile'].write(s + '\n')  
    
    # chrom to do separate
    #int names
    chromNames = []
    chromNames.append('chrX')  # put X first, so start with X, chr1 -- longest chroms start first

    for i in range(1,39):
        c = 'chr' + str(i)
        chromNames.append(c)
    
    s = 'have %i chroms to do individually' % len(chromNames)
    print(s,flush=True)
    myData['logFile'].write(s + '\n')  
    
    toDoTogether = []
    for c in contigNames:
        if c not in chromNames:
            toDoTogether.append(c)
    
    s = 'have %i chroms to do together' % len(toDoTogether)
    print(s,flush=True)
    myData['logFile'].write(s + '\n')  
    
    # make intervals files
    myData['bqsrIntervalsDir'] = myData['workingBaseDir'] + 'bqsr-intervals'
    
    if os.path.isdir(myData['bqsrIntervalsDir']) is True:
        s = '%s exists!' % (myData['bqsrIntervalsDir'])
        print(s,flush=True)
        myData['logFile'].write(s + '\n')
    else:
       cmd = 'mkdir %s' % myData['bqsrIntervalsDir']
       myData['logFile'].write(cmd + '\n')
       runCMD(cmd)    
    myData['bqsrIntervalsDir'] += '/'
           
    myData['bqsrIntervalsFiles'] = []    
    for c in chromNames:
        fn = myData['bqsrIntervalsDir'] + c + '.list'
        outFile = open(fn,'w')
        outFile.write('%s\n' % c)
        outFile.close()
        myData['bqsrIntervalsFiles'].append(fn)
        
    # then add the rest
    fn = myData['bqsrIntervalsDir'] + 'other' + '.list'
    outFile = open(fn,'w')
    for c in toDoTogether:
        outFile.write('%s\n' % c)
    outFile.close()
        
    myData['bqsrIntervalsFiles'].insert(0,fn) # so that starts with the other, this is large want more run time
        
    s = 'have setup %i bqsrIntervalsFiles files' % len(myData['bqsrIntervalsFiles'])    
    print(s,flush=True)
    myData['logFile'].write(s + '\n')  
    myData['logFile'].flush()  
############################################################################# 

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

# if log file exists, then there is partial processing so not sure we want to redo and overwrite
# safe to just quite and letter user deal with it

#if os.path.isfile(myData['logFileName']) is True:
#    print('ERROR!!!')
#    print('%s exists.  Do you really want to rerun this pipeline?' % myData['logFileName'] )
#    print('ERROR!!!')
#    sys.exit()

myData['logFile'] = open(myData['logFileName'],'w')

# add initial infoto log
init_log(myData)

# make sure programs are availble
check_prog_paths(myData)
check_dir_space(myData)
run_bwa_mem2(myData,False)
run_mdspark(myData,False)

run_bqsr(myData,False)



# clean up and get elapsed time!
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



