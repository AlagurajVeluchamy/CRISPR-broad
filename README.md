# CRISPR-broad
**CRISPR-broad**
**Citation**
Alaguraj Veluchamy, Kaian Teles and Wolfgang Fischle. CRISPR-broad: combined design of multitargeting gRNAs and broad, multiplex target finding (Sci Rep. 2023 Nov 12;13(1):19717.).**

**https://pubmed.ncbi.nlm.nih.gov/37953351/**


Version: 1.1.0

Date: 20-12-2022

Authors: Alaguraj Veluchamy alaguraj.v@gmail.com, Wolfgang Fischle


The package is developed for the design of gRNA for the targeted epigenome modifications on a broader region.

License: GPL-3

**Depends:**
bwa-0.7.11, 
Python 3.9.1, 
biopython 1.78,
pandas 1.2.0
pyranges-0.0.95

**Installing bwa:**
git clone https://github.com/lh3/bwa.git; 
cd bwa; make

**Installing dependancies through pip:**
pip install biopython;
pip install pandas;
pip install pyranges


**Installing CRISPR-broad:**

git clone https://github.com/AlagurajVeluchamy/CRISPR-broad.git; cd CRISPR-broad

**Description:**

CRISPR-broad is a standalone tool that  enables user to scan genome for regions that has high frequency of gRNA with user-supplied variation. This open source tool is robust and efficient in finding gRNAs with scores and ranks potential target region. 

**Input:**

CRISPR-broad minimally only requires a genome file in fasta format, PAM sequences in string format "NGG".

**QUICK RUN:**

If you are planning to split genome into 40 parts in slurm based cluster: srun --time=4:00:00 --mem=100G --nodes=1 -c 40  --pty bash -l 

1. python crisprbroad.py genomesplit -d /Users/CRISPR -f Inputgenome.fa -t 40 -l 23  -g 50 -p GG

2. python crisprbroad.py createindex -f Inputgenome.fa

3. python crisprbroad.py maptogenome -d /Users/CRISPR -f Inputgenome.fa -t 40 -nm 5 -nx 10000 -m 2 -g 50 -l 230

4. python crisprbroad.py filterhits -d /Users/CRISPR  -t 40 -nm 5 -nx 10000

6. python crisprbroad.py findwindow -d /Users/CRISPR  -t 40  -p GG -w 10000 -l 23 -nm 5 -nw 5


**Usage: To list all modules**

usage: python crisprbroad.py [-h] 

genomesplit,  createindex,  maptogenome,  filterhits,  findwindow, findmultiwindow, multisgrna

**Module 1: Split genome into candidate gRNA sequences**

python crisprbroad.py genomesplit -h

    usage: crisprbroad.py genomesplit [-h] 
                                          -h, --help            show this help message and exit
                                          -f GENOMESPLITFASTA, --genome_fasta GENOMESPLITFASTA
                                                                Genome sequence in FASTA format
                                          -q INPUTBED           --query_bedfile query in bed format
                                                                optional input for finding sgRNA in a given region
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

python crisprbroad.py createindex -h

    usage: python crisprbroad.py createindex [-h] 
                                                 -h, --help            show this help message and exit
                                                 -f GENOMESPLITFASTA, --genome_fasta GENOMESPLITFASTA
                                                                Genome sequence in FASTA format

**Module 3: Map gRNA to genome**

python crisprbroad.py maptogenome -h
    
    usage: crisprbroad.py maptogenome [-h]
                                          -h, --help            show this help message and exit
                                          -f GENOMESPLITFASTA, --genome_fasta GENOMESPLITFASTA
                                                                Genome sequence in FASTA format
                                          -d WORKINGDIRECTORY, --working_directory WORKINGDIRECTORY
                                                                Complete path of output directory
                                          -m MISMATCH, --num_mismatch MISMATCH
                                                                Maximum number of mismatches for Bowtie2 alignment
                                          -t THREADS, --num_threads THREADS
                                                                Launch t number of threads in parallel
                                          -nm MINHITS, --get_minhits MINHITS
                                                                Minimum number of hits in a window (filter)
                                          -nx MAXHITS, --get_maxmum_total_hits MAXHITS
                                                                Maximum total number of hits
                                          -g GC, --get_gc GC    gc content in integer
                                          -l CANDIDATERNALENGTH, --get_candidaternalength CANDIDATERNALENGTH
                                                                Candidate gRNA length
                                       
**Module 4:Filter hits for on-target and off-target**

python crisprbroad.py filterhits -h

    usage: crisprbroad.py filterhits [-h]
                                          -h, --help            show this help message and exit
                                          -d WORKINGDIRECTORY, --working_directory WORKINGDIRECTORY
                                                                Complete path of output directory
                                          -t THREADS, --num_threads THREADS
                                                                Launch t number of threads in parallel
                                          -nm MINHITS, --get_minhits MINHITS
                                                                Minimum number of hits in a window (filter)
                                          -nx MAXHITS, --get_maxmum_total_hits MAXHITS

**Module 5: Scoring windows and ranking gRNA**

python crisprbroad.py findwindow -h

    usage: crisprbroad.py findwindow [-h] 
                                          -h, --help            show this help message and exit
                                          -as AVOIDSCOREOFF, --avoid_scoreoff AVOIDSCOREOFF
                                                                Avoid score calculation for off-targets (0 or 1)
                                          -d WORKINGDIRECTORY, --working_directory WORKINGDIRECTORY
                                                                Complete path of output directory
                                          -p PAMSEQUENCE, --pam_sequence PAMSEQUENCE
                                                                PAM sequence string eg: NGG
                                          -t THREADS, --num_threads THREADS
                                                                Launch t number of threads in parallel
                                          -nm MINHITS, --get_minhits MINHITS
                                                                Minimum number of hits in a window (filter)
                                          -w WINDOW, --get_window WINDOW
                                                                Window size in bp
                                          -l CANDIDATERNALENGTH, --get_candidaternalength CANDIDATERNALENGTH
                                                                Candidate gRNA length
                                          -nw WINDOWNUMBERS, --number of window WINDOWNUMBERS
                                                               number of windows for targeting

**Module 6: gRNA that targets multiple windows**

python crisprbroad.py findmultiwindow -h

    usage: crisprbroad.py findmultiwindow [-h] 
                                          -h, --help            show this help message and exit
                                          -f GENOMESPLITFASTA, --genome_fasta GENOMESPLITFASTA
                                                                Genome sequence in FASTA format
                                          -d WORKINGDIRECTORY, --working_directory WORKINGDIRECTORY
                                                                Complete path of output directory
                                          -p PAMSEQUENCE, --pam_sequence PAMSEQUENCE
                                                                PAM sequence string eg: NGG
                                          -t THREADS, --num_threads THREADS
                                                                Launch t number of threads in parallel
                                          -nm MINHITS, --get_minhits MINHITS
                                                                Minimum number of hits in a window (filter)
                                          -nx MAXHITS, --get_maxmum_total_hits MAXHITS
                                                                Maximum total number of hits
                                          -ws WINDOW, --get_window WINDOW
                                                                Window size in bp
                                          -sl SLIDINGWINDOWSIZE, --get_slidingwindowlength SLIDINGWINDOWSIZE
                                                                sliding window in bp (For consecutive windows: same as window size)
                                          -nw NUMBEROFWINDOW, --get_multiwindow WINDOW
                                                                number of target windows
                                          -q INPUTBED           --query_bedfile query in bed format
                                                                optional input for finding sgRNA in a given region
                                          -l CANDIDATERNALENGTH, --get_candidaternalength CANDIDATERNALENGTH
                                                                Candidate gRNA length

**Module 7: Finding more than one gRNA that targets multiple windows**

python crisprbroad.py multisgrna -h

    usage: crisprbroad.py multisgrna [-h] 
                                          -h, --help            show this help message and exit
                                          
                                          -d WORKINGDIRECTORY, --working_directory WORKINGDIRECTORY
                                                                Complete path of output directory
                                          -cb CRISPRBROADRESULTFILE, --crisprbroadresultfilt CRISPRBROADRESULTFILE
                                                                File obtained as result from module5
                                          -mg SGRNANUMBERS, --numberofgRNA NUMBEROFSGRNA
                                                                Number of sgRNA in target window
  


**RESULTS:**

**OUTPUT description **

Results ranked list of gRNA and windows with score and dispersion are in tab-limited file: CRISPR_broad_results.xls**



| Column       |    Column Name        | Description |
| ------------- |:-------------:|------------- |
|1| Chromosome      | Chromosome name of best window for gRNA selection | 
|2| Start     | Start position of best window for gRNA selection | 
|3| End     | End position of best window for gRNA selection | 
|4| crnaid     | Assigned Id for sgRNA (forward/reverse strand; sgRNA number; Start position) | 
|5| Sequence     | sgRNA sequence of user defined length    | 
|6| stringforgroup |  Chr-st-end pattern |
|7| Hits within Best window      | Higher number of hits => Better sgRNA      |
|8| score     | Score assigned (0-2) for sgRNA (Higher the score better is the candidate |
|9| stdev     | Dispersion/spread from the center of the hits      |
|10| Finalscore     | Sum of score of all hits sgRNA (Higher the score better is the candidate |
|11| Total hits |Number of sgRNA hits in  whole genome including off-targets  | 
|12| offtargethits     | Number of sgRNA hits in  outside Best target window based on score |       |



**CIGAR Description:**

Example format: ["('I', '+', '12072082', '23M', '0')"]
The above is repeated (equal to the number of sgRNA hits)

| CIGAR string       |    Desction        | 
| ------------- |:-------------:|
| I | chromosome name |
| +       |    strand        | 
| 12072082 | start |
| 23M | 23 nucleotide match i,e alignment of sgRNA (23nt length) |
| 0 | Edit distance (Unused) |
