#!/usr/bin/env python3

'''
CRISP v1.0: Collapse/Retain Investigation of SNP Phylogenetic 
by
Tomi Jacobs
'''
#pip install ipynb-py-convert
import sys
import os
import shutil
import subprocess


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


# STEP 2: BUILD HISAT2 INDEX
#https://daehwankimlab.github.io/hisat2/manual/
def build_index(trans, outdir):
    ''' this builds a hisat index from each transcriptome'''
    index_dir = os.path.join(outdir, "index")
    index_prefix = os.path.join(index_dir, "transcriptome")
    if os.path.exists(index_prefix + ".1.ht2"):
      print("index exists, skipping...")
      return index_prefix  
    os.makedirs(index_dir, exist_ok=True)
    cmd = "hisat2-build " + trans + " " + index_prefix
    run(cmd, "hisat2-build")
    return index_prefix

# MAIN: this runs only when script is executed directly
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("USE: python crisp.py transcriptome.fa r1.fq r2.fq output.fa")
        sys.exit(1)

    trans  = sys.argv[1]
    r1     = sys.argv[2]
    r2     = sys.argv[3]
    output = sys.argv[4]

    for f in [trans, r1, r2]:
        if not os.path.exists(f):
            sys.exit("file not found: " + f)

    if trans.endswith(".gz"):
        sys.exit("transcriptome must be zipped: " + trans)

    workdir = output.replace(".fa", "_workdir").replace(".fasta", "_workdir")
    os.makedirs(workdir, exist_ok=True)

    check_dependencies()
    index_prefix = build_index(trans, workdir)