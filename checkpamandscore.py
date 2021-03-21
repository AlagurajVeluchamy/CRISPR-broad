import os
import re
import string
from Bio import SeqIO
import subprocess
from subprocess import call, Popen, PIPE
import glob
from Bio.Seq import Seq
##from Bio.Alphabet import generic_dna
import numpy
import pandas as pd
import pyranges as pr
import sys
import csv
csv.field_size_limit(sys.maxsize)



def checkpam(argu, chr,start,secondaryhitpattern,strand, getgenome, querysequence, secondaryhitmismatch):
    length = score = ggfound = 0
    forwardseq = revcomplementseq = complementseq = "NONE"
    for record in getgenome:
        recordid = str(record.id)
        #if (strand == 0 or strand <= 256) and recordid == chr:
        if (strand == "+" ) and recordid == chr:
            forwardseq = str(record.seq)[start-1:start+22]
            if forwardseq.endswith(argu.pamsequence):
                ggfound = ggfound + 1
            elif argu.pamsequence in forwardseq[-4:]:
                ggfound = ggfound + 1
        elif (strand == "-") and recordid == chr:
            complementseq = str(record.seq)[start-1:start + 22]
            newseq = Seq(complementseq)
            revcomplementseq = str(newseq.reverse_complement())
            if revcomplementseq.endswith(argu.pamsequence):
                ggfound = ggfound + 1
            elif argu.pamsequence in revcomplementseq[-4:]:
                ggfound = ggfound + 1
    if ggfound > 0:
        #print ("inside checkpam.. before getscore")
        (length, score) = getscore(secondaryhitpattern, secondaryhitmismatch)
    return(length, score, ggfound,forwardseq,revcomplementseq)

def getscore(firstlineunused, pattern1,secondaryhitmismatch):
    pattern = str(pattern1)
    #print(pattern)
    #print("patternprinted...")
    #print (secondaryhitmismatch)
    #secondaryhitmismatch = int(secondaryhitmismatch1)
    #print (pattern)
    #print (secondaryhitmismatch)
    matchpat = re.findall("(\d+)M", pattern)
    #mismatchpat = re.findall("(\d+)X", pattern)
    insertionpat = re.findall("(\d+)I", pattern)
    deletionpat = re.findall("(\d+)D", pattern)
    substitutionpat = re.findall("(\d+)D", pattern)
    #print("inside checkpam.. inside getscore")
    if len(matchpat) > 0:
        summatchpat = sum(map(int, matchpat))
    else:
        summatchpat = 0
    if len(substitutionpat) > 0:
        substitutionpat = sum(map(int, substitutionpat))
    else:
        substitutionpat = 0
    if len(insertionpat) > 0:
        suminsertionpat = sum(map(int, insertionpat))
    else:
        suminsertionpat = 0
    length = summatchpat + substitutionpat + suminsertionpat
    #score = round(float(summatchpat - secondaryhitmismatch - suminsertionpat)/float(length), 3)
    score = (summatchpat - secondaryhitmismatch - suminsertionpat) / length
    ##return(length, score)
    #print(score)
    return (score)


def calculatewrite(argu, fileoutmap, gc_content, counthits, countoffhits, inwindscore, outwindscore, queryname, querychr, arbitraryid, querystart, offtargetlist, secondaryhitschr, bestwindow, newlist, querysequence, getstddev):
    import string
    suminwindscore = sum(map(float, inwindscore))
    sumoutwindscore = sum(map(float, outwindscore))
    finalscore = ((float(counthits) * (suminwindscore/float(counthits)))/float(counthits+countoffhits)) - ((float(countoffhits) * (sumoutwindscore/float(countoffhits)))/float(counthits+countoffhits))
    ### Below sequence translate take more time ###
    if counthits > argu.minhits:
        if "candidaternars" in queryname:
            ## For python 3.0 it changed to str.maketrans
            #complements = string.maketrans('ATGC', 'TACG')
            #complements = str.maketrans('ATGC', 'TACG')
            #print (querysequence)
            #querysequence = querysequence.translate(complements)[::-1]
            #print (querysequence)
            ############ Above sequence translate take more time ###
            fileoutmap.write(str(queryname) + "_" + str(querychr) + "_" + str(arbitraryid) + "_" + str(querystart) + "\t" + str(querysequence) + "\t" + str(gc_content)+ "\t" + str(finalscore) + "\t" + str(getstddev) + "\t" + str(suminwindscore) + "\t" + str(sumoutwindscore) +  "\t" + str(counthits) + "\t" + str(len(offtargetlist)) + "\t" + str(querychr) + ":" + str(min(bestwindow)) + "-" + str(max(bestwindow)) + "\t" + str(newlist) + "\n")
        elif "candidaternafs" in queryname:
            fileoutmap.write(str(queryname) + "_" + str(querychr) + "_" + str(arbitraryid) + "_" + str(querystart) + "\t" + str(querysequence) + "\t" + str(gc_content) + "\t" + str(finalscore) + "\t" + str(getstddev) + "\t" + str(suminwindscore) + "\t" + str(sumoutwindscore) + "\t" + str(counthits) + "\t" + str(len(offtargetlist)) + "\t" + str(querychr) + ":" + str(min(bestwindow)) + "-" + str(max(bestwindow)) + "\t" + str(newlist) + "\n")
    elif counthits < 1:
        print ("Pattern Match error")

def calculatestdev(bestwindow):
    if len(bestwindow) > 4:
        getstddev = numpy.std(bestwindow, axis=0)
    else:
        getstddev = "NULL"
    return getstddev

def getallresults(argu):
    dataresult = pd.read_csv("Crispr_broad_multihits.xls", names=['stringforgroup', 'hitsinwindow', 'score', 'stdev', 'finalscore','crnaid', 'Sequence', 'Totalhits','offtargetshits', 'Chromosome', 'Start', 'End'], sep="\t", header=None)
    print ("Pooling the results...")
    column_names = ['Chromosome', 'Start', 'End',  'crnaid', 'Sequence', 'stringforgroup', 'hitsinwindow', 'score', 'stdev', 'finalscore', 'Totalhits', 'offtargetshits' ]
    dataresult = dataresult.reindex(columns=column_names)
    dataresult1 = dataresult.sort_values(['offtargetshits', 'crnaid'], ascending=[True, True])
    dataresult1.to_csv("Crispr_broad_multihits_sorted.xls", sep="\t", index=False, header=True)


def checkpandapamquerywindow (rowdatchunk, argu, splitinputfile):
    genomesplitoriginal = "Candidate_crisprnatabbed.txt"
    datachunkallcrna = pd.read_csv(genomesplitoriginal, header=None, names=['name', 'Str', 'Chromosome', 'crnanumber', 'Startggposition'], sep="\t", dtype=str, engine='c', index_col=False, low_memory=False, na_filter=True)
    eachrowdataframe = pd.DataFrame(rowdatchunk['colCIGAR'], columns=['flips'])
    eachrowdataframe['flips'] = eachrowdataframe['flips'].str.replace('XA:Z:', '', regex=False)
    eachrowdataframe['flips'] = eachrowdataframe['flips'].str.replace('+', '+,', regex=False)
    eachrowdataframe['flips'] = eachrowdataframe['flips'].str.replace('-', '-,', regex=False)
    eachrowdataframe = eachrowdataframe['flips'].str.replace('XA:Z:', '', regex=False)
    eachrowdataframe1 = pd.DataFrame(eachrowdataframe.str.split(',').tolist(), columns=['Chromosome', 'Str', 'Start', 'mismatch', 'unusedstrings'])

    eachrowdataframe2 = eachrowdataframe1.copy()
    eachrowdataframe2.loc[eachrowdataframe1.index.max()+1] = [rowdatchunk['chr'],rowdatchunk['strand'], rowdatchunk['firsthitstart'], '23M','0']
    eachrowdataframe2 = eachrowdataframe2.dropna()

    eachrowdataframe2['score'] = eachrowdataframe2.apply(lambda x: getscore(x, x['mismatch'], pd.to_numeric(x['unusedstrings']) ), axis=1)
    eachrowdataframe2.loc[:,'taggingcount'] = 1
    #########################

    eachrowdataframe3 = eachrowdataframe2[['Chromosome', 'Start', 'score','taggingcount']]
    eachrowdataframe3.loc[:,'End'] = eachrowdataframe3['Start'].astype(int) + 27
    eachrowdataframe3 = eachrowdataframe3.dropna()

    datachunkallcrna1 = pd.DataFrame()
    datachunkallcrna1['Chromosome'] = datachunkallcrna['Chromosome']
    datachunkallcrna1['Start'] = datachunkallcrna['Startggposition']
    datachunkallcrna1.loc[:,'Start'] = datachunkallcrna1['Start'].astype(int) - 4
    datachunkallcrna1.loc[:,'End'] = datachunkallcrna['Startggposition'].astype(int) + 4

    overlappedggdataframe = (pr.PyRanges(eachrowdataframe3).overlap(pr.PyRanges(datachunkallcrna1))).df

    if (len(overlappedggdataframe) >= int(argu.minhits) & len(overlappedggdataframe) <= int(argu.maxhits)):
        #allchrmax = overlappedggdataframe[['Chromosome', 'Start']].groupby('Chromosome').Start.agg('max').reset_index()
        #allchrmin = overlappedggdataframe[['Chromosome', 'Start']].groupby('Chromosome').Start.agg('min')
        #allchrnew = pd.merge(allchrmin, allchrmax, on='Chromosome')
        #allchrnew.columns = ['Chromosome', 'Start', 'End']
        ################### Window from the given ########
        filemakewindowsinwindow = pd.read_csv(argu.inputbed, sep="\t", names=['Chromosome', 'Start', 'End'])
        filemakewindowsinwindow1 =pr.PyRanges(filemakewindowsinwindow)
        ###################
        getscorewindow = (filemakewindowsinwindow1.join(pr.PyRanges(overlappedggdataframe))).df
        suminwindscore = getscorewindow['score'].sum()
        counthits = getscorewindow['score'].shape[0]
        getscorewindow['stringforgroup'] = getscorewindow[['Chromosome', 'Start', 'End']].astype(str).agg('-'.join, axis=1)
        #getscorewindow1 = getscorewindow.groupby(['stringforgroup'], as_index=False).agg({'score': "sum"})
        #print (getscorewindow)
        getscorewindowsum = getscorewindow.groupby(['stringforgroup'], as_index=False)['score'].sum()
        #print(getscorewindow2)
        getscorewindowcount = getscorewindow.groupby(['stringforgroup'], as_index=False)['taggingcount'].count()
        getscorewindow2 = pd.merge(getscorewindowcount, getscorewindowsum, on='stringforgroup')
        #print (getscorewindow2)
        hitslistforstdevdf = getscorewindow.groupby('stringforgroup')['Start_b'].apply(list).reset_index(name='starthistaslist')
        #print(hitslistforstdevdf)
        getscorewindow2["stdev"] = hitslistforstdevdf.apply(lambda x: calculatestdev( x['starthistaslist']), axis=1)
        #print (getscorewindow2)

        #getscorewindow1[['Chromosome', 'Start', 'End']] = getscorewindow1['stringforgroup'].str.split("-", expand=True)
        #print(getscorewindow2)
        ##########
        getscorenonwindow = (pr.PyRanges(overlappedggdataframe).join(filemakewindowsinwindow1, how="left")).df
        sumoutwindscore = getscorenonwindow[getscorenonwindow['Start_b'] < 0]['score'].sum()
        countoffhits = getscorenonwindow[getscorenonwindow['Start_b'] < 0].shape[0]
        #############

        #suminwindscore = sum(map(float, inwindscore))
        #sumoutwindscore = sum(map(float, outwindscore))
        finalscore = ((float(counthits) * (suminwindscore / float(counthits))) / float(counthits + countoffhits)) - (
                    (float(countoffhits) * (sumoutwindscore / float(countoffhits))) / float(counthits + countoffhits))
        ###########
        #print (finalscore)
        #getstddev = calculatestdev(bestwindow)
        getscorewindow2.loc[:,"finalscore"] =  finalscore

        getscorewindow2["crnaid"] = rowdatchunk['crnaid']
        getscorewindow2["sequence"] = rowdatchunk['sequence']
        getscorewindow2.loc[:,"Totalhits"] = len(overlappedggdataframe)
        #print(getscorewindow2)
        ####################
        getscorewindow3 = getscorewindow2.sort_values(['taggingcount'], ascending=[False]).head(argu.windownumbers)
        # print(data1['Totalhits'].iloc[0])
        # print(data1)
        getscorewindow3['offtargetshits'] = pd.to_numeric(getscorewindow3['Totalhits']).iloc[0] - getscorewindow3['taggingcount'].sum()
        # data1['off-targets-hits'] = data1.apply(lambda row: TotalHitsforthatRNA - inwindowtarget, axis=1)
        # print ("I am writing")
        getscorewindow3[['Chromosome', 'Start', 'End']] = getscorewindow3['stringforgroup'].str.split("-", expand=True)
        #print (getscorewindow3)
        fileoutstd = open(os.path.join(argu.workingdirectory, "Crispr_broad_multihits.xls"), 'a')
        getscorewindow3.to_csv(fileoutstd, sep="\t", index=False, header=False)

def checkpandapam(rowdatchunk, argu, splitinputfile):
    genomesplitoriginal = "Candidate_crisprnatabbed.txt"
    datachunkallcrna = pd.read_csv(genomesplitoriginal, header=None, names=['name', 'Str', 'Chromosome', 'crnanumber', 'Startggposition'], sep="\t", dtype=str, engine='c', index_col=False, low_memory=False, na_filter=True)
    eachrowdataframe = pd.DataFrame(rowdatchunk['colCIGAR'], columns=['flips'])
    eachrowdataframe['flips'] = eachrowdataframe['flips'].str.replace('XA:Z:', '', regex=False)
    eachrowdataframe['flips'] = eachrowdataframe['flips'].str.replace('+', '+,', regex=False)
    eachrowdataframe['flips'] = eachrowdataframe['flips'].str.replace('-', '-,', regex=False)
    eachrowdataframe = eachrowdataframe['flips'].str.replace('XA:Z:', '', regex=False)
    eachrowdataframe1 = pd.DataFrame(eachrowdataframe.str.split(',').tolist(), columns=['Chromosome', 'Str', 'Start', 'mismatch', 'unusedstrings'])

    eachrowdataframe2 = eachrowdataframe1.copy()
    eachrowdataframe2.loc[eachrowdataframe1.index.max()+1] = [rowdatchunk['chr'],rowdatchunk['strand'], rowdatchunk['firsthitstart'], '23M','0']
    eachrowdataframe2 = eachrowdataframe2.dropna()

    eachrowdataframe2['score'] = eachrowdataframe2.apply(lambda x: getscore(x, x['mismatch'], pd.to_numeric(x['unusedstrings']) ), axis=1)
    eachrowdataframe2.loc[:,'taggingcount'] = 1
    #########################

    eachrowdataframe3 = eachrowdataframe2[['Chromosome', 'Start', 'score','taggingcount']]
    #eachrowdataframe3['End'] = eachrowdataframe3['Start'].astype(int) + 27
    eachrowdataframe3.loc[:,'End'] = eachrowdataframe3['Start'].astype(int) + 27
    eachrowdataframe3 = eachrowdataframe3.dropna()
    #print (eachrowdataframe3)

    datachunkallcrna1 = pd.DataFrame()
    datachunkallcrna1['Chromosome'] = datachunkallcrna['Chromosome']
    datachunkallcrna1['Start'] = datachunkallcrna['Startggposition']
    datachunkallcrna1.loc[:,'Start'] = datachunkallcrna1['Start'].astype(int) - 4
    datachunkallcrna1.loc[:,'End'] = datachunkallcrna['Startggposition'].astype(int) + 4

    overlappedggdataframe = (pr.PyRanges(eachrowdataframe3).overlap(pr.PyRanges(datachunkallcrna1))).df

    if (len(overlappedggdataframe) >= int(argu.minhits) & len(overlappedggdataframe) <= int(argu.maxhits)):
        allchrmax = overlappedggdataframe[['Chromosome', 'Start']].groupby('Chromosome').Start.agg('max').reset_index()
        allchrmin = overlappedggdataframe[['Chromosome', 'Start']].groupby('Chromosome').Start.agg('min')

        allchrnew = pd.merge(allchrmin, allchrmax, on='Chromosome')
        allchrnew.columns = ['Chromosome', 'Start', 'End']
        filemakewindowsinwindow1 = pr.PyRanges(allchrnew).window(argu.windowsize)
        #filemakewindowsinwindow1 = pr.PyRanges(filemakewindowsinwindow)

        getscorewindow = (filemakewindowsinwindow1.join(pr.PyRanges(overlappedggdataframe))).df
        suminwindscore = getscorewindow['score'].sum()
        counthits = getscorewindow['score'].shape[0]
        getscorewindow['stringforgroup'] = getscorewindow[['Chromosome', 'Start', 'End']].astype(str).agg('-'.join, axis=1)
        #getscorewindow1 = getscorewindow.groupby(['stringforgroup'], as_index=False).agg({'score': "sum"})
        #print (getscorewindow)
        getscorewindowsum = getscorewindow.groupby(['stringforgroup'], as_index=False)['score'].sum()
        #print(getscorewindow2)
        getscorewindowcount = getscorewindow.groupby(['stringforgroup'], as_index=False)['taggingcount'].count()
        getscorewindow2 = pd.merge(getscorewindowcount, getscorewindowsum, on='stringforgroup')
        #print (getscorewindow2)
        hitslistforstdevdf = getscorewindow.groupby('stringforgroup')['Start_b'].apply(list).reset_index(name='starthistaslist')
        #print(hitslistforstdevdf)
        getscorewindow2["stdev"] = hitslistforstdevdf.apply(lambda x: calculatestdev( x['starthistaslist']), axis=1)
        print (getscorewindow2)

        #getscorewindow1[['Chromosome', 'Start', 'End']] = getscorewindow1['stringforgroup'].str.split("-", expand=True)
        #print(getscorewindow2)
        ##########
        getscorenonwindow = (pr.PyRanges(overlappedggdataframe).join(filemakewindowsinwindow1, how="left")).df
        sumoutwindscore = getscorenonwindow[getscorenonwindow['Start_b'] < 0]['score'].sum()
        countoffhits = getscorenonwindow[getscorenonwindow['Start_b'] < 0].shape[0]
        #############

        #suminwindscore = sum(map(float, inwindscore))
        #sumoutwindscore = sum(map(float, outwindscore))
        finalscore = ((float(counthits) * (suminwindscore / float(counthits))) / float(counthits + countoffhits)) - (
                    (float(countoffhits) * (sumoutwindscore / float(countoffhits))) / float(counthits + countoffhits))
        ###########
        #print (finalscore)
        #getstddev = calculatestdev(bestwindow)
        getscorewindow2.loc[:,"finalscore"] =  finalscore

        getscorewindow2["crnaid"] = rowdatchunk['crnaid']
        getscorewindow2["sequence"] = rowdatchunk['sequence']
        getscorewindow2.loc[:,"Totalhits"] = len(overlappedggdataframe)
        #print(getscorewindow2)
        ####################
        getscorewindow3 = getscorewindow2.sort_values(['taggingcount'], ascending=[False]).head(argu.windownumbers)
        # print(data1['Totalhits'].iloc[0])
        # print(data1)
        getscorewindow3['offtargetshits'] = pd.to_numeric(getscorewindow3['Totalhits']).iloc[0] - getscorewindow3['taggingcount'].sum()
        # data1['off-targets-hits'] = data1.apply(lambda row: TotalHitsforthatRNA - inwindowtarget, axis=1)
        # print ("I am writing")
        getscorewindow3[['Chromosome', 'Start', 'End']] = getscorewindow3['stringforgroup'].str.split("-", expand=True)
        #print (getscorewindow3)
        fileoutstd = open(os.path.join(argu.workingdirectory, "Crispr_broad_multihits.xls"), 'a')
        getscorewindow3.to_csv(fileoutstd, sep="\t", index=False, header=False)

def overlapeachchrommultiwindow(argu, splitinputfile):
    print ("Processing:",splitinputfile)
    for datachunk in pd.read_csv(splitinputfile, header=None, names =['crnaid', 'strand', 'chr', 'firsthitstart', 'sequence', 'colX0', 'colX1', 'colCIGAR','firsthits', 'Totalhits'], sep="\t", comment='@', dtype=str, engine='c', index_col=False, low_memory=False, na_filter=True, chunksize=200000, iterator=True):
        datachunk["colCIGAR"] = datachunk["colCIGAR"].str.split(";", expand=False)
        datachunk = datachunk[datachunk["colCIGAR"].str.len() > int(argu.minhits)]
        #datachunk["colCIGAR"] = datachunk["colCIGAR"].str.split(";", expand=False)
        pd.set_option('display.max_rows', 10)
        #print (datachunk)
        #print (datachunk["colCIGAR"])
        #datachunk.apply(lambda x: print (x['crnaid']), axis =1)
        if argu.inputbed is None:
            datachunk.apply(lambda x: checkpandapam(x, argu, splitinputfile), axis =1 )
        else:
            datachunk.apply(lambda x: checkpandapamquerywindow(x, argu, splitinputfile), axis=1)

def multisgrna(argu):
    #multisgrnafile = "crispr_broad_multisgrna.xls"
    multisgrnafile= pr.PyRanges(pd.read_csv(argu.crisprbroadresfile, sep="\t"))
    getscorewindow = (multisgrnafile.join(multisgrnafile)).df
    #overlappedggdataframe = (pr.PyRanges(eachrowdataframe2).overlap(pr.PyRanges(datachunkallcrna1))).df
    getscorewindow.to_csv(os.path.join(argu.workingdirectory, "crispr_broad_multisgrna.xls"), sep="\t", index=False, header=True)
    #filemultisgrnain = open(os.path.join(argu.workingdirectory, "crispr_broad_multisgrna.xls"), 'w+')
    #p2 = subprocess.Popen(['bedtools', 'intersect', '-a', str(argu.crisprbroadresfile), '-b', str(argu.crisprbroadresfile), '-wao'], stdout=filemultisgrnain)
    #p2.communicate()
