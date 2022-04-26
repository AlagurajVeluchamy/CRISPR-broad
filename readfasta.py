from Bio import SeqIO
import re
from subprocess import call
import glob
import os
import math
import pyranges as pr

def readgenomefastafile(argu, inputfastafile):
    """ Read genome file in fasta format and split into multiple files. """
    print ("Genome file:", inputfastafile)
    print ("Chosen PAM sequence :", argu.pamsequence)
    fileoutcandidates = os.path.join(argu.workingdirectory, "Candidate_crisprna.txt")
    fileoutcandidates1 = open(fileoutcandidates, 'w')
    count = 0
    for record in SeqIO.parse(inputfastafile, "fasta"):
        print ("Finding PAM pattern in the chromosome :", str(record.id))
        forwardseq = str(record.seq)
        start = re.finditer(argu.pamsequence,str(record.seq))
        complementseq = str(record.seq.complement())
        complementedstart =  re.finditer(argu.pamsequence,complementseq)
        for i in start:
            pos = i.start()
            candidaterna = forwardseq[pos-(argu.candidaternalength-2):pos+2]
            if len(candidaterna) > 0:
                gc_content = (float((candidaterna.count('G') + candidaterna.count('C'))) / len(candidaterna)) * 100
                if gc_content > argu.gc:
                    fileoutcandidates1.write('>cfs_' +str(record.id)+'_'+ str(count) + "_"+ str(pos+2)+'\n')
                    fileoutcandidates1.write(candidaterna + '\n')
                    count += 1
        for j in complementedstart:
            poscomp = j.start()
            candidaternareverserna = complementseq[poscomp:poscomp+argu.candidaternalength][::-1]
            if len(candidaternareverserna) > 0:
                gc_content = (float((candidaternareverserna.count('G') + candidaternareverserna.count('C'))) / len(candidaternareverserna)) * 100
                if gc_content > argu.gc:
                    fileoutcandidates1.write('>crs_' +str(record.id)+'_'+ str(count) + "_"+ str(poscomp+1)+'\n')
                    fileoutcandidates1.write(candidaternareverserna + '\n')
                    count += 1
    return count

def createtabbedpamfile(argu):
    """ Read genome file in fasta format and split into multiple files. """
    print ("Create PAM positions in tabbed format:", argu.genomesplitfasta)
    #print ("Chosen PAM sequence :", argu.pamsequence)
    fileoutcandidatestabbed = os.path.join(argu.workingdirectory, "Candidate_crisprnatabbed.txt2col")
    fileoutcandidates1tabbed = open(fileoutcandidatestabbed, 'w')
    count = 0
    for record in SeqIO.parse(argu.genomesplitfasta, "fasta"):
        forwardseq = str(record.seq)
        start = re.finditer(argu.pamsequence,str(record.seq))
        complementseq = str(record.seq.complement())
        complementedstart =  re.finditer(argu.pamsequence,complementseq)
        for i in start:
            pos = i.start()
            candidaterna = forwardseq[pos-(argu.candidaternalength-2):pos+2]
            if len(candidaterna) > 0:
                fileoutcandidates1tabbed.write(str(record.id) + "\t" +  str(pos + 2) + '\n')
                count += 1
        for j in complementedstart:
            poscomp = j.start()
            candidaternareverserna = complementseq[poscomp:poscomp+argu.candidaternalength][::-1]
            if len(candidaternareverserna) > 0:
                fileoutcandidates1tabbed.write(str(record.id) + "\t" + str(poscomp+1) + '\n')
                count += 1
    return count


def readbedtofasta(argu):
    inputbedpr = pr.read_bed(str(os.path.join(argu.workingdirectory, argu.inputbed)))
    print (str(os.path.join(argu.workingdirectory,argu.genomesplitfasta)))
    inputbedseq = pr.get_fasta(inputbedpr, str(os.path.join(argu.workingdirectory,argu.genomesplitfasta)))
    convertquerybedtofasta = os.path.join(argu.workingdirectory, argu.inputbed+".fa")
    fileoutinputbedseq = open(convertquerybedtofasta, 'w')
    for eachbedlineseq, eachbedlinechr, eachbedlinestart, eachbedlineend in zip(inputbedseq, inputbedpr.df.Chromosome, inputbedpr.df.Start, inputbedpr.df.End ):
        #print (eachbedline[1].Chromosome)
        fileoutinputbedseq.write('>' + str(eachbedlinechr) + "KKK" + str(eachbedlinestart) + "KKK" + str(eachbedlineend) + '\n' + str(eachbedlineseq) + '\n')
    return convertquerybedtofasta
    #print (inputbedseq)
    #readgenomefastafile(argu.)
    #print (argu.genomesplitfasta)
    #print (str(os.path.join(argu.workingdirectory, argu.genomesplitfasta)))
    #pr.read_bed(filesinput)
    #pr.read_bed(path, nrows=5)

def splitinputfileformulti(argu,countinputlines):
    print ("Total Candidate gRNA in the genome:", countinputlines)
    print ("Splitting the genome for multiprocessing into: ", int(argu.threads), "files")
    print ("Splitting the genome for multiprocessing: Number of sequences per process = ", str(math.ceil((countinputlines/int(argu.threads))/2.0) * 2))
    #call(['split', '-l', str(int(2*(countinputlines/int(argu.threads)))), '--additional-suffix=.crisprfastq', str(os.path.join(argu.workingdirectory,'Candidate_crisprna.txt')), str(os.path.join(argu.workingdirectory,'Candidate_crisprna.split'))])
    call(['split', '-l', str(int((math.ceil((countinputlines/int(argu.threads))/2.0) * 4)+2)), str(os.path.join(argu.workingdirectory,'Candidate_crisprna.txt')),  str(os.path.join(argu.workingdirectory,'Candidate_crisprna.split'))])
    splitinputfiles = []
    for file in glob.glob(str(os.path.join(argu.workingdirectory,"Candidate_crisprna.split*"))):
        splitinputfiles.append(file)
    return splitinputfiles

def getminhits(argu):
    getminhitspattern = ""
    for minofhit in range(argu.minhits-1):
        getminhitspattern = getminhitspattern + str(".*;")
    return getminhitspattern


def getlistfiles(argu):
    print (argu.workingdirectory)
    filesinput = glob.glob(str(os.path.join(argu.workingdirectory, "Candidate_crisprna.split")) + "*")
    filteredfilename = []
    for fileinput in filesinput:
        #print (fileinput)
        if not ".txt"  in fileinput and not ".sam" in fileinput:
            filteredfilename.append(fileinput)
            print (fileinput)
    return filteredfilename

def getsamlistfiles(argu):
    filesaminput = glob.glob(str(os.path.join(argu.workingdirectory, "Candidate_crisprna.split")) + "*" + "sam")
    return filesaminput

def getsamfilteredfiles(argu):
    filteredsamfilename = glob.glob(str(os.path.join(argu.workingdirectory, "Candidate_crisprna.split")) + "*" + ".sam_filtered.txt")
    return filteredsamfilename
