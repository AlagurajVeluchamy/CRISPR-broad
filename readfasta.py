import os
from argumentparse import *
from Bio import SeqIO
import re

def readgenomefastafile(argu):
    """ Read bed/bam files from a directory. """
    print argu.genomesplitfasta
    print argu.pamsequence
    fileoutcandidates = os.path.join(argu.outputdir, "Candidate_crisprna.txt")
    fileoutcandidates1 = open(fileoutcandidates, 'w')
    for record in SeqIO.parse(argu.genomesplitfasta, "fasta"):
        print(record.id)
        print len(str(record.seq))
        forwardseq = str(record.seq)
        start = re.finditer("GG",str(record.seq))
        #Total number of matches start1 = re.findall("GG", str(record.seq))
        complementseq = str(record.seq.complement())
        complementedstart =  re.finditer("GG",complementseq)
        count =0
        for i in start:
            count+=1
            pos = i.start()
            candidaterna = forwardseq[pos-(argu.candidaternalength-2):pos+2]
            fileoutcandidates1.write('>candidaternafs_' +str(record.id)+'_'+ str(count) + "_"+ str(pos+2)+'\n')
            fileoutcandidates1.write(candidaterna + '\n')
        for j in complementedstart:
            count+=1
            poscomp = j.start()
            candidaternareverserna = complementseq[poscomp:poscomp+argu.candidaternalength][::-1]
            fileoutcandidates1.write('>candidaternars_' +str(record.id)+'_'+ str(count) + "_"+ str(poscomp+1)+'\n')
            fileoutcandidates1.write(candidaternareverserna + '\n')

            #print candidaterna
        #print(record.seq[start - 23:start+2])
    return file

