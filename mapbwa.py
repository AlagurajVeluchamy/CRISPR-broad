from subprocess import call
import os
import pandas as pd
import re
import string
from itertools import islice

def mapbwafastatogenome(argu):
    ##Indexing the genome
    call(['bwa', 'index', argu.genomesplitfasta])
    call(['bwa', 'aln', '-t', argu.threads, '-n', argu.mismatch,  '-f', str(os.path.join(argu.outputdir,'Candidate_crisprna.aln')), argu.genomesplitfasta, str(os.path.join(argu.outputdir,'Candidate_crisprna.txt'))])
    call(['bwa', 'samse', '-n', '1000', '-f', str(os.path.join(argu.outputdir,'Candidate_crisprna.sam')), argu.genomesplitfasta, str(os.path.join(argu.outputdir,'Candidate_crisprna.aln')), str(os.path.join(argu.outputdir,'Candidate_crisprna.txt'))])

def filtersam(argu):
    colnames = ['crnaid', 'waste1', 'chr', 'hitstart', 'waste2', 'cigar']
    data = pd.read_csv(os.path.join(argu.outputdir, 'Candidate_crisprna.sam'), sep="\t", names=colnames, usecols=['crnaid', 'waste1', 'chr', 'hitstart', 'waste2', 'cigar', 'waste3', 'waste4', 'waste5', 'sequence'], skiprows=8, dtype=str)
    data['merged'] = data['chr'].astype(str) + '_' + data['hitstart'].astype(str) + '_' + data['cigar'].astype(str)
    data1 = data[['crnaid', 'sequence', 'merged']].groupby('crnaid').agg({'merged': lambda x: ','.join(x)}).reset_index()
    datafiltered = data1[data1['merged'].str.contains(",", na=False)]
    datafiltered.to_csv(os.path.join(arguoutputdir, 'Candidate_crisprna_filtered.txt'), sep="\t", index=False, header=False)
    ###working data = pd.read_csv(os.path.join(argu.outputdir, 'Candidate_crisprna.sam'), sep="\t", names=(range(20)))
    ###working datafiltered = data[data[19].str.contains("XA:Z", na=False)]
    ###working datafiltered.to_csv(os.path.join(argu.outputdir, 'Candidate_crisprna_filtered.txt'), sep="\t", index=False, header=False)


def overlapeachchromosome(argu):
    """
    :type argu: directory
    """
    fileinmap = open(os.path.join(argu.outputdir, 'Candidate_crisprna_filtered.txt'), 'r')
    fileoutmap = open(os.path.join(argu.outputdir, 'Crisprna_result.xls'), 'w')
    fileoutmap.write("AQueryname\tQuery Sequence\tGC %\tNumber of hits in Best window\tNumber of hits outside Best window\tBest window\tHits within Best window\tHits outside Best window\n")
    complements = string.maketrans('ATGC', 'TACG')
    for line in fileinmap:
        line1 = re.search("(\w*)_(\w*)_(\w*)_(\d*)\t(\w*)\t(\w*)\t(\w*)\t(\w.*)\t(\w*)\t(.*?)\t(\w.*)\t(\w.*)\t(\w*)\t(.*?).*XA:Z:(\w*.*)", str(line))
        if line1:
            #print line1.group(1)
            newlist = []
            offtargetlist = []
            (queryname, querychr, arbitraryid, querystart, firsthitchr, firsthitstart, querysequence) = line1.group(1,2,3,4,6,7,13)
            gc_content = float((querysequence.count('G') + querysequence.count('C'))) / len(querysequence) * 100
            line2 = re.findall("(\w*.*?),[+-](\d*.*?),.*?;", str(line1.group(15)))
            counthits = 0
            bestwindow = []
            for eachtuple in line2:
                if querychr == eachtuple[0] and ((int(querystart) - argu.window) <= int(eachtuple[1]) <= (int(querystart) + argu.window)):
                    bestwindow.append(int(eachtuple[1]))
                    counthits = counthits + 1
                    newlist.append(str(eachtuple))
                else:
                    offtargetlist.append(str(eachtuple))
            if counthits > argu.minhits:
                ### Below sequence translate take more time ###
                if "candidaternars" in queryname:
                    querysequence = querysequence.translate(complements)[::-1]
                ############ Above sequence translate take more time ###
                fileoutmap.write(str(queryname) + "_" + str(querychr) + "_" + str(arbitraryid) + "_" + str(querystart) + "\t" + str(querysequence) + "\t" + str(gc_content) + "\t" + str(counthits) + "\t"+str(len(offtargetlist)) + "\t" + str(eachtuple[0])+":"+str(min(bestwindow)) + "-" + str(max(bestwindow)) + "\t" + str(newlist) + "\n")
        else:
            print "Pattern Match error"
#     if ($seq_s){print FASTA "${id}_S\n$seq_s\n";}
#     if ($seq_a){
#     $seq_a=~tr / atucgACGUT / TAAGCTGCAA /;
#     $seq_a=reverse($seq_a);
#     print FASTA "${id}_A\n$seq_a\n";
