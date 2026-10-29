#!/usr/bin/env python3

'''
Project 2
CRISP v1.0: Collapse/Retain Investigation of SNP Phylogenetic with Bowtie2
by
Tomi Jacobs, Jordan Kressman & Joe Walker

What we run on terminal: python collab_crisp.py test/ChloroSubSet_hetSites.fasta test/ChloroSubset_1.fq test/ChloroSubset_2.fq ChloroSubset_hetSites_masked.fa
python collab_crisp.py test/ChloroSubSet_hetSites.fasta test/ChloroSubSet_1.fq test/ChloroSubSet_2.fq ChloroSubSet_hetSites_masked.fa
For introducing heterozygosity: python introduce_het.py test/ChloroSubSet.fasta test/ChloroSubSet_hetSites.fasta test/het_sites_log.txt
sudo python collabcrisp_bowtie.py ChloroSubSet_hetSites.fasta test/ChloroSubSet_1.fq test/ChloroSubSet_2.fq ChloroSubset_hetSites_Bowtie2Masked.fa

Shawn files trial

When fastq is in data
python collabcrisp_bowtie.py \
    ~/data/projectTwo/tomi_all_fasta/Wrightia_religiosa.transdecoder.cds.fasta \
    ~/data/projectTwo/jordan_all_fastq/Wrightia_religiosa_1.fq \
    ~/data/projectTwo/jordan_all_fastq/Wrightia_religiosa_2.fq \
    ~/data/projectTwo/CRISP_bowtie/44fasta_maskedwithbowtie/Wrightia_religiosa_masked.fasta

When fastq is in scratch
python collabcrisp_bowtie.py ~/data/projectTwo/tomi_all_fasta/Tabernaemontana_elegans.transdecoder.cds.fasta /scratch/jordirectory/Tabernaemontana_elegans_1.fq /scratch/jordirectory/Tabernaemontana_elegans_2.fq ~/data/projectTwo/CRISP_bowtie/44fasta_maskedwithbowtie/Tabernaemontana_elegans_masked.fasta
'''

import sys
import os
import shutil
import subprocess

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python crisp.py transcriptome.fa fastq1.fq fastq2.fq output.fa")
        sys.exit(1)
    transcriptome = sys.argv[1]
    r1            = sys.argv[2]
    r2            = sys.argv[3]
    output        = sys.argv[4]

    # CHECK DEPENDENCIES
    #["bowtie2", "bowtie2-build", "samtools", "bcftools"]
    samtoolspath = shutil.which("samtools")
    if samtoolspath is None:
        print("samtools not detected")
        sys.exit(1)
    else:
        print("samtools detected: success")

    bowtie2buildpath = shutil.which("bowtie2-build")

    if bowtie2buildpath is None:
        print("bowtie2-build not detected")
        sys.exit(1)
    else:
        print("bowtie2-build detected: success")

    bowtie2path = shutil.which("bowtie2")
    if bowtie2path is None:
        print("bowtie2 not detected")
        sys.exit(1)
    else:
        print("bowtie2 detected: success")

    bcftoolspath = shutil.which("bcftools")
    if bcftoolspath is None:
        print("bcftools not detected")
        sys.exit(1)
    else:
        print("bcftools detected: success")

#MAIN

#STEP 1 BUILD BOWTIE INDEX 2
#Uisng Burrows-Wheeler Transform, hisat2 rearranges the content of our transcriptome file (assembled by TJ,hetsites by Jordan)
#and rearranges all its content into a special sorted structure(). 
# #-c path to the fasta files
#bowtie2-build [options]* <reference_in> <bt2_base>
#if os.path.exists("bowtie_index.1.bt2"):
    #print("step 1 done, skipping...")
#else:
print("\nStep 1: building bowtie2 index..")
subprocess.run(["bowtie2-build", "--threads", "10", transcriptome, "bowtie_index"])
print("bowtie2 index built, success")
#hisat2-build [options]* <reference_in> <ht2_index>  


#STEP 2 ALIGN READS TO TRANSCRIPTOME
#hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]
#hisat2 usage follows this format
#align and saves sam file
#its no spliced cause its transcriptome and not genome
#-x: index prefix
#-1 & -2: paired end reads r1 & r2
#-S: output SAM file
#bowtie2 -x cp_index -1 reads_R1.fastq -2 reads_R2.fastq -S aligned_output
#if os.path.exists("output.bam"):
#   print("step 2 done, skipping...")
#else:
print("\nStep 2: Aligning reads to transcriptome")
subprocess.run(["bowtie2", "-p", "10", "-x","bowtie_index","-1",r1,"-2",r2,"-S", "output.sam"])
subprocess.run(["samtools", "sort", "-@", "10", "-o", "output.bam", "output.sam"])
    #convert sam(large plain text file) to bam(compressed smaller file)
    # #-o output bam file
    # # subprocess.run(["samtools", "sort", "-o", "output.bam", "output.sam"])
    # #index the BAM file, creates output output.bam.bai file which allows tools jump positions quickly 
subprocess.run(["samtools", "index", "output.bam"])
    #we can ddelete the sam file since we no longer nead it
subprocess.run(["rm", "output.sam"])
print("alignment complete: success")


#STEP 3: ASSESS ALL COVERAGE
#if os.path.exists("coverage.depth"):
#   print("step 3 done, skipping...")
#else:
print("\nStep 3: assessing coverage...")
subprocess.run("samtools depth -a output.bam > coverage.depth", shell=True)
print("coverage assessment complete: success")
#samtools depth helps counts how many reads cover each position
#-a this helps us read all coverage including zeros.Output: transcript, position & depth


#STEP 4: CALL AND FILTER SNP/VARIANTS --BCF tools 
#All SNPs are genetic variants but not all variants are SNPs
#bcfools mpileup stacks all reads at every position.. what nucleotide each read shows at every position
#variat sites are where bfttools identifies positions where reads disagree 
#mv outputting only variant sites, -OC uncompressed VCF output format
#if os.path.exists("filtered_variants.vcf"):
#   print("step 4 done, skipping...")
#else:
print("\nStep 4: calling variants...")
subprocess.run("bcftools mpileup --threads 10 -f " + transcriptome + " output.bam | bcftools call -mv -Ov -o raw_variants.vcf", shell=True)
subprocess.run("bcftools filter -i 'QUAL>=30 && DP>=10 && MQ>=20' raw_variants.vcf > filtered_variants.vcf", shell=True)
    
        #lets count the number of variants/snps that survive the filter so we know just for refernce sake
count = sum(1 for line in open("filtered_variants.vcf") if not line.startswith ("#"))
print("variants passing filter: " + str(count))
print("variant calling complete: success")

#to keep only high confidence variant/hetrozygous sites
#Phred scores
#QUAL (Variant Quality score)>= 30: 99.9% its not a sequencing error but a real varinat 
#DP (Read Depth)>= 10: minimum 10 reads to cover this position
#MQ (Mapping Quality)>= 20: filters ambiguiosly mapped reads
   
#STEP 5: MASK HETEROZYGOUS SITES ***
#if os.path.exists(output):
#    print("step 5 done, skipping...")
#else:
print("\nStep 5: masking heterozygous sites...")
vcf = open("filtered_variants.vcf","r")
het_sites = open(transcriptome,"r")
newfile = open(output,"w")

vcf_dictionary = {}
for line in vcf:
    if "#" in line:
        pass
    else:
        linesplit = line.split('\t')
        if len(linesplit[3]) != 1 or len(linesplit[4]) != 1:
            pass
        else:
            if (linesplit[0]) in vcf_dictionary:
                vcf_dictionary[(linesplit[0])].append(linesplit[1])
            else:
                vcf_dictionary[(linesplit[0])] = [linesplit[1]]

line = ''
chromy = []
stupidvariable = 0
linetitle = ''

for line in het_sites:
    if stupidvariable == 0:
        stupidvariable += 1
        linetitle = (((line.split(" "))[0]).strip(">"))
        newfile.write(line)
    else:
        if ">" in line:
            try:
                for i in vcf_dictionary[linetitle]:
                    chromy[int(i)-1] = "N"
            except KeyError:
                pass
            sillystring = "".join(chromy)
            for k in range(0, len(sillystring), 60):
                newfile.write(sillystring[k:k+60] + "\n")
            newfile.write(line)
            linetitle = (((line.split())[0]).strip(">"))
            chromy = []
        else:
            chromy += (list(line.strip("\n")))

try:
    for i in vcf_dictionary[linetitle]:
        chromy[int(i)-1] = "N"
except KeyError:
    pass
sillystring = "".join(chromy)
for k in range(0, len(sillystring), 60):
    newfile.write(sillystring[k:k+60] + "\n")

print("masking complete: success")


blockart = r"""
    __  ____   ____ _____ ____          __   ___   ___ ___  ____  _        ___ ______    ___  ___   
   /  ]|    \ |    / ___/|    \        /  ] /   \ |   |   ||    \| |      /  _]      |  /  _]|   \  
  /  / |  D  ) |  (   \_ |  o  )      /  / |     || _   _ ||  o  ) |     /  [_|      | /  [_ |    \ 
 /  /  |    /  |  |\__  ||   _/      /  /  |  O  ||  \_/  ||   _/| |___ |    _]_|  |_||    _]|  D  |
/   \_ |    \  |  |/  \ ||  |       /   \_ |     ||   |   ||  |  |     ||   [_  |  |  |   [_ |     |
\     ||  .  \ |  |\    ||  |       \     ||     ||   |   ||  |  |     ||     | |  |  |     ||     |
 \____||__|\_||____|\___||__|        \____| \___/ |___|___||__|  |_____||_____| |__|  |_____||_____|
                                                                                                    
"""
print(blockart)

print("output written to file:",output)
