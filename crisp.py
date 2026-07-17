#!/usr/bin/env python3

'''
CRISP v1.0: Collapse/Retain Investigation of SNP Phylogenetic 
by
Tomi Jacobs

This pipeline would:
1. Map reads to transcriptomes using HISAT2
2. Calculate coverage and identify low-coverage regions
3. Mask ambiguous regions (low coverage) as 'N'
4. Call SNPs using bcftools
5. Generate collapsed and retained sequences???
'''

import sys
import os
import csv
import shutil
import subprocess
from pathlib import Path


# CHECK DEPENDENCIES
def check_dependencies():
    #a list of all external tools needed for crisp to run. mafft, iq tree were taken out.
    #limiting this script to the masking
    tools = ["hisat2", "hisat2-build", "samtools", "bcftools"]
    missing = []
    for t in tools:
        if shutil.which(t) is None:
            missing.append(t)

    #missing = [t for t in tools if not shutil.which(t)] #this shows what tool is missing

    if missing:
        #gives instruction on what tool should be installed amongst the list
        sys.exit(f"Missing: {', '.join(missing)}\nPlease install these tools first")
    print("All dependencies found")

check_dependencies()


# HELPER: RUN SHELL COMMAND
def run(cmd, desc):
    '''Executes shell command with progress tracking'''
    #prints description of the step currently running- for tracking progress
    print("running: " + desc)
    #execute shell command sing subprocess
    result = subprocess.run(cmd, shell=True)
    #check if command failed
    if result.returncode != 0:
        #this describes the actual step that failed
        print("failed: " + desc)
        #this prints the command that failed 
        print("command was: " + cmd)
        sys.exit() #end
    print("done: " + desc) #success


# STEP 1: VALIDATE SAMPLESHEET (TEST SAMPLE)
def read_samplesheet(samplesheet):
    '''Reads and validates the samplesheet'''
    samples = []
    if not os.path.exists(samplesheet):
        sys.exit("samplesheet not found: " + samplesheet)

    f = open(samplesheet, "r")
    reader = csv.DictReader(f)
    for col in ["sample", "r1", "r2", "transcriptome"]:
        if col not in reader.fieldnames:
            sys.exit("missing column: " + col)
    for row in reader:
        sample = row["sample"].strip()
        r1     = row["r1"].strip()
        r2     = row["r2"].strip()
        trans  = row["transcriptome"].strip()
        for filepath in [r1, r2, trans]:
            if not os.path.exists(filepath):
                sys.exit("file not found: " + filepath)
        if trans.endswith(".gz"):
            sys.exit("transcriptome must be unzipped: " + trans)
        samples.append([sample, r1, r2, trans])
    f.close()
    if len(samples) == 0:
        sys.exit("no samples found in samplesheet")
    print("samplesheet ok - " + str(len(samples)) + " samples found")
    return samples


# STEP 2: BUILD HISAT2 INDEX
#https://daehwankimlab.github.io/hisat2/manual/
def build_index(samples, outdir):
    ''' this builds a hisat index from each transcriptone'''

    print ("\n ---step 2: building hisat2 indices---")

    for sample, r1, r2, trans in samples:
        index_dir = os.path.join(outdir, sample, "index")
        index_prefix = os.path.join(index_dir, sample)

        if os.path.exists(index_prefix + ".1.ht2"):
            print("index exists, skip: " + sample)
            continue

        os.makedirs(index_dir, exist_ok=True)

        cmd = "hisat2-build " + trans + " " + index_prefix
        run(cmd, "hisat-build for " + sample)


# MAIN: this runs only when script is executed directly
if __name__ == "__main__":

    #check arguments
    if len(sys.argv) != 3:
        print("USE: python TJ_crisp.py samplesheet.csv outdir/")
        sys.exit(1)

    samplesheet = sys.argv[1]
    outdir      = sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    #run pipeline
    check_dependencies()
    samples = read_samplesheet(samplesheet)