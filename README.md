# EpiCRISPR-TargetFinder
**EpiCRISPR-TargetFinder**

Version: 1.1.0

Date: 20-05-2020

Authors: Alaguraj Veluchamy alaguraj.v@gmail.com, Wolfgang Fischle


The package is developed for the design of gRNA for the targeted epigenome modifications on a broader region.

License: GPL-2

**Depends:**
bwa-0.7.11, 
Python 3.9.1, 
biopython 1.78,
pandas 1.2.0

**Installing bwa:**
git clone https://github.com/lh3/bwa.git; 
cd bwa; make

**Installing dependancies through pip:**
pip install biopython;
pip install pandas;


**Installing EpiCRISPR-TargetFinder:**

git clone https://github.com/AlagurajVeluchamy/EpiCRISPR-TargetFinder.git; cd EpiCRISPR-TargetFinder

**Description:**

EpiCRISPR-targetfinder is a standalone tool that  enables user to scan genome for regions that has high frequency of gRNA with user-supplied variation. This open source tool is robust and efficient in finding gRNAs with scores and ranks potential target region. 

**Input:**

EpiCRISPR-targetfinder minimally only requires a genome file in fasta format, PAM sequences in string format "NGG".

**QUICK RUN:**

python epicrisprtarget.py genomesplit -d /Users/EpiCRISPR -f Inputgenome.fa -g 50 -p GG -t 1 -l 23

python epicrisprtarget.py createindex -f Inputgenome.fa

python epicrisprtarget.py maptogenome -d /Users/EpiCRISPR -f Inputgenome.fa -m 1 -n 3 -k 1 -g 50 -l 23 -t 1

python epicrisprtarget.py filterhits -d /Users/EpiCRISPR -n 3 -t 1

python epicrisprtarget.py findwindow -d /Users/EpiCRISPR -f Inputgenome.fa -p GG -t 1 -l 23 -w 1000 -n 3


**Usage: To list all modules**

usage: python epicrisprtarget.py [-h] 

genomesplit,

createindex,

maptogenome,

filterhits,

findwindow

**Module 1: Split genome into candidate gRNA sequences**

python epicrisprtarget.py genomesplit -h

    usage: epicrisprtarget.py genomesplit [-h] 
                                          -h, --help            show this help message and exit
                                          -f GENOMESPLITFASTA, --genome_fasta GENOMESPLITFASTA
                                                                Genome sequence in FASTA format
                                          -d WORKINGDIRECTORY, --working_directory WORKINGDIRECTORY
                                                                Complete path of output directory
                                          -p PAMSEQUENCE, --pam_sequence PAMSEQUENCE
                                                                PAM sequence string eg: GG for NGG PAM
                                          -l CANDIDATERNALENGTH, --get_candidaternalength CANDIDATERNALENGTH
                                                                Candidate gRNA length
                                          -g GC, --get_gc GC    gc content in integer
                                          -t THREADS, --num_threads THREADS
                                                                Launch t number of threads in parallel

**Module 2: Create index**

python epicrisprtarget.py createindex -h

    usage: python epicrisprtarget.py createindex [-h] 
                                                 -h, --help            show this help message and exit
                                                 -f GENOMESPLITFASTA, --genome_fasta GENOMESPLITFASTA
                                                                Genome sequence in FASTA format

**Module 3: Map gRNA to genome**

python epicrisprtarget.py maptogenome -h
    
    usage: epicrisprtarget.py maptogenome [-h]
                                          -h, --help            show this help message and exit
                                          -f GENOMESPLITFASTA, --genome_fasta GENOMESPLITFASTA
                                                                Genome sequence in FASTA format
                                          -d WORKINGDIRECTORY, --working_directory WORKINGDIRECTORY
                                                                Complete path of output directory
                                          -m MISMATCH, --num_mismatch MISMATCH
                                                                Maximum number of mismatches for Bowtie2 alignment
                                          -t THREADS, --num_threads THREADS
                                                                Launch t number of threads in parallel
                                          -n MINHITS, --get_minhits MINHITS
                                                                Minimum number of hits in a window (filter)
                                          -k MAXHITS, --get_maxmum_total_hits MAXHITS
                                                                Maximum total number of hits
                                          -g GC, --get_gc GC    gc content in integer
                                          -l CANDIDATERNALENGTH, --get_candidaternalength CANDIDATERNALENGTH
                                                                Candidate gRNA length
                                       
**Module 4:Filter hits for on-target and off-target**
python epicrisprtarget.py filterhits -h

    usage: epicrisprtarget.py filterhits [-h]
                                          -h, --help            show this help message and exit
                                          -d WORKINGDIRECTORY, --working_directory WORKINGDIRECTORY
                                                                Complete path of output directory
                                          -t THREADS, --num_threads THREADS
                                                                Launch t number of threads in parallel
                                          -n MINHITS, --get_minhits MINHITS
                                                                Minimum number of hits in a window (filter)

**Module 5: Scoring windows and ranking gRNA**
python epicrisprtarget.py findwindow -h

    usage: epicrisprtarget.py findwindow [-h] 
                                          -h, --help            show this help message and exit
                                          -f GENOMESPLITFASTA, --genome_fasta GENOMESPLITFASTA
                                                                Genome sequence in FASTA format
                                          -d WORKINGDIRECTORY, --working_directory WORKINGDIRECTORY
                                                                Complete path of output directory
                                          -p PAMSEQUENCE, --pam_sequence PAMSEQUENCE
                                                                PAM sequence string eg: NGG
                                          -t THREADS, --num_threads THREADS
                                                                Launch t number of threads in parallel
                                          -n MINHITS, --get_minhits MINHITS
                                                                Minimum number of hits in a window (filter)
                                          -w WINDOW, --get_window WINDOW
                                                                Window size in bp
                                          -l CANDIDATERNALENGTH, --get_candidaternalength CANDIDATERNALENGTH
                                                                Candidate gRNA length


**RESULTS:
Results ranked list of gRNA and windows with score and dispersion are in tab-limited file: EpiCRISPR_results.xls**
