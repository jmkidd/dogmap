import sys
import glob
import argparse


###############################################################################

# SETUP

parser = argparse.ArgumentParser(description='come stats')

parser.add_argument('--dir', type=str,help='dir containing *stats.txt output',required=True)
parser.add_argument('--out', type=str,help='name of file for writing stats table',required=True)
args = parser.parse_args()

#####################################################################



fn = args.out

flList = glob.glob(args.dir + '/*stats.txt')
print('have %i for stats' % len(flList))

flList.sort()


outFile = open(fn,'w')
nl = []
inFile = open(flList[0],'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    nl.append(line[0])
inFile.close()
nl = '\t'.join(nl) + '\n'
outFile.write(nl)



for f in flList:
    nl = []
    inFile = open(f,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        nl.append(line[1])
    inFile.close()
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)


print('DONE!')
print('stats written to',fn)
outFile.close()