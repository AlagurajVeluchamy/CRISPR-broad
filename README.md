# EpiCRISPR-TargetFinder
**EpiCRISPR-TargetFinder**

Version: 1.1.0

Date: 20-05-2020

Authors: Alaguraj Veluchamy alaguraj.v@gmail.com, Wolfgang Fischle


The package is developed for the design of gRNA for the targeted epigenome modifications on a broader region.

License: GPL-2

Depends: bwa-0.7.17-r1188, python3, biopython, pandas, sort
Installing bwa:
git clone https://github.com/lh3/bwa.git
cd bwa; make

Installing dependancies through pip:
pip install biopython
pip install pandas


**Installation:**
EpiCRISPR-TargetFinder depends on python3, bwa. 

**Description:**
EpiCRISPR-targetfinder is a standalone tool that  enables user to scan genome for regions that has high frequency of gRNA with user-supplied variation. This open source tool is robust and efficient in finding gRNAs with scores and ranks potential target region. 

**Input:**
EpiCRISPR-targetfinder minimally only requires a genome file in fasta format, PAM sequences in string format "NGG".

**Usage: To list all modules**
usage: python epicrisprtarget.py [-h] 

**Module 1: Split genome into candidate gRNA sequences**

python epicrisprtarget.py genomesplit -h

    usage: epicrisprtarget.py genomesplit -d GENOMESPLITFASTA
                                          -p PAMSEQUENCE
                                          -o OUTPUTDIR
                                          -m MISMATCH
                                          -t THREADS
                                          -s GETSEQUENCES
                                          -n MINHITS
                                          -w WINDOW
                                          -l CANDIDATERNALENGTH

**Module 2: Create index**

python epicrisprtarget.py createindex -h

    usage: python epicrisprtarget.py createindex [-h] -f GENOMESPLITFASTA

**Module 3: Map gRNA to genome**

python epicrisprtarget.py maptogenome -h
    
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
python epicrisprtarget.py filterhits -h

    usage: epicrisprtarget.py filterhits [-h]
                                          -d WORKINGDIRECTORY
                                          -t THREADS
                                          -n MINHITS

**Module 5: Scoring windows and ranking gRNA**
python epicrisprtarget.py findwindow -h

    usage: epicrisprtarget.py findwindow [-h] 
                                          -f GENOMESPLITFASTA 
                                          -d WORKINGDIRECTORY
                                          -p PAMSEQUENCE 
                                          -t THREADS
                                          -n MINHITS 
                                          -w WINDOW
                                          -l CANDIDATERNALENGTH

**RESULTS:
Results ranked list of gRNA and windows with score and dispersion are in tab-limited file: EpiCRISPR_results.xls**

Arguments list:
    
      -f GENOMESPLITFASTA, --genome_fasta GENOMESPLITFASTA
                            Genome sequence in FASTA format
      -c CHR, --chromosome CHR
                            chromosome name
      -d WORKINGDIRECTORY, --working_directory WORKINGDIRECTORY
                            Complete path of output directory
      -p PAMSEQUENCE, --pam_sequence PAMSEQUENCE
                            PAM sequence string eg: NGG
      -l CANDIDATERNALENGTH, --get_candidaternalength CANDIDATERNALENGTH
                            Candidate gRNA length
      -g GC, --get_gc GC    gc content in integer
      -t THREADS, --num_threads THREADS
                            Launch t number of threads in parallel
      -n MINHITS, --get_minhits MINHITS
                            Minimum number of hits in a window (filter)
      -w WINDOW, --get_window WINDOW
                            Window size in bp

