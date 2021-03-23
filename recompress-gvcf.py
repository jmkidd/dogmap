# recompress-gvcf.py

# program to recompress a gvcf file to save space


import sys
import subprocess
import os
import argparse
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
    for p in ['bgzip','tabix']:
        if shutil.which(p) is None:
            s = p + ' not found in path! please fix (module load?)'
            print(s, flush=True)
            sys.exit()
        else:
            print('%s\t%s' % (p,shutil.which(p)),flush=True)

###############################################################################

# SETUP

parser = argparse.ArgumentParser(description='recompress gvcf, leaves indicator file when done')

parser.add_argument('--gvcf', type=str,help='sample gvcf file',required=True)

args = parser.parse_args()

#####################################################################

myData = {} # dictionary for keeping and passing information
myData['gvcf'] = args.gvcf
myData['gvcfCompleteToken'] = myData['gvcf'] + '.recompressed'

if os.path.isfile(myData['gvcfCompleteToken']) is True:
    print('token exists, exiting\n%s\n' % myData['gvcfCompleteToken'] )
    sys.exit()

# make sure programs are availble
check_prog_paths(myData)

# original size
print('\nStarting\n')
myData['origSize'] = os.path.getsize(myData['gvcf'])
myData['origSize'] = myData['origSize'] / (1024**3) # to GB

print('File name: %s' % myData['gvcf'])
print('Original size (Gb): %.3f' % (myData['origSize']))




# compress
tmpNewGvcf = myData['gvcf'] + '.tmp.gz'

cmd = 'zcat %s | bgzip -@ 6 -l 9  -c > %s' % (myData['gvcf'],tmpNewGvcf)
print(cmd,flush=True)
runCMD(cmd)

# rename
cmd = 'mv %s %s ' % (tmpNewGvcf,myData['gvcf'])
print(cmd,flush=True)
runCMD(cmd)

# index
cmd = 'tabix -p vcf -f %s ' % myData['gvcf']
print(cmd,flush=True)
runCMD(cmd)

# get new size
myData['newSize'] = os.path.getsize(myData['gvcf'])
myData['newSize'] = myData['newSize'] / (1024**3) # to GB

print('File name: %s' % myData['gvcf'])
print('New size (Gb): %.3f' % (myData['newSize']))

f = myData['newSize'] / myData['origSize']
print('New is %.2f %% of the original ' % (f * 100.0) )

cmd = 'touch %s' % myData['gvcfCompleteToken']
print(cmd,flush=True)
runCMD(cmd)

print('\nComplete!\n',flush=True)
