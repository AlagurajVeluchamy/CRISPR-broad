# EpiCRISPR-TargetFinder
**EpiCRISPR-TargetFinder**

Version: 1.1.0

Date: 20-05-2020

Authors: Alaguraj Veluchamy alaguraj.v@gmail.com, Wolfgang Fischle


The package is developed for the design of gRNA for the targeted epigenome modifications on a broader region.

License: GPL-2

Depends: bowtie-2, pandas

**Installation:**
EpiCRISPR-TargetFinder depends on python3, bowtie. 

**Description:**
EpiCRISPR-targetfinder is a standalone tool that  enables user to scan genome for regions that has high frequency of gRNA with user-supplied variation. This open source tool is robust and efficient in finding gRNAs with scores and ranks potential target region. 

**Input:**
EpiCRISPR-targetfinder minimally only requires a genome file in fasta format, PAM sequences in string format "NGG".

**Usage:**
usage: epicrisprtarget.py [-h] 

**Module 1: Split genome into candidate gRNA sequences**

python epicrisprtarget.py genomesplit -h

usage: epicrisprtarget.py genomesplit [-d GENOMESPLITFASTA]
                                      [-p PAMSEQUENCE] [-o OUTPUTDIR]                                      
                                      [-m MISMATCH] [-t THREADS]                                      
                                      [-s GETSEQUENCES] [-n MINHITS]                                   
                                      [-w WINDOW] [-l CANDIDATERNALENGTH]

**Module 2: Create index**
usage: epicrisprtarget.py createindex [-h] -f GENOMESPLITFASTA

**Module 3: Map gRNA to genome**
usage: epicrisprtarget.py maptogenome [-h] 
                                       -f GENOMESPLITFASTA
                                       -d WORKINGDIRECTORY 
                                       -m MISMATCH 
                                       -t THREADS
                                       -s GETSEQUENCES 
                                       -n MINHITS 
                                       -g GC content in %
                                       -l CANDIDATERNALENGTH
                                       
**Module 4:Filter hits for on-target and off-target**
usage: epicrisprtarget.py filterhits [-h]
                                      -d WORKINGDIRECTORY 
                                      -t THREADS 
                                      -n MINHITS

**Module 5: Scoring windows and ranking gRNA**
usage: epicrisprtarget.py findwindow [-h] 
                                      -f GENOMESPLITFASTA 
                                      -d WORKINGDIRECTORY 
                                      -p PAMSEQUENCE 
                                      -t THREADS 
                                      -n MINHITS 
                                      -w WINDOW 
                                      -l CANDIDATERNALENGTH

Arguments:
    -h, --help            show this help message and exit

    -d GENOMESPLITFASTA, --genome_fasta GENOMESPLITFASTA
                          Read bed/bam files from a directory.

    -p PAMSEQUENCE, --pam_sequence PAMSEQUENCE
                          Read bed/bam files from a directory.

    -o OUTPUTDIR, --output_directory OUTPUTDIR
                          Read bed/bam files from a directory.

    -m MISMATCH, --num_mismatch MISMATCH
                          Read bed/bam files from a directory.

    -t THREADS, --num_threads THREADS
                          Read bed/bam files from a directory.

    -s GETSEQUENCES, --get_sequence GETSEQUENCES
                          Read bed/bam files from a directory.

    -n MINHITS, --get_minhits MINHITS
                          Read bed/bam files from a directory.

    -w WINDOW, --get_window WINDOW
                          Read bed/bam files from a directory.

    -l CANDIDATERNALENGTH, --get_candidaternalength CANDIDATERNALENGTH
                          Read bed/bam files from a directory.

