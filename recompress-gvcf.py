# recompress-gvcf.py

# program to recompress a gvcf file to save space


import sys
import subprocess
import os
import argparse
import shutil

import recompressgvcf




if __name__ == '__main__':
    # SETUP
    parser = argparse.ArgumentParser(description='recompress gvcf, leaves indicator file when done')
    parser.add_argument('--gvcf', type=str,help='sample gvcf file',required=True)
    args = parser.parse_args()    
    myData = {} # dictionary for keeping and passing information
    myData['gvcf'] = args.gvcf
    
    recompressgvcf.run_recompress_gvcf(myData)


