import os
import re
import string
from Bio import SeqIO
import subprocess
from subprocess import call, Popen, PIPE
import glob
from Bio.Seq import Seq
##from Bio.Alphabet import generic_dna
import numpy as np
import numpy
import pandas as pd
import pyranges as pr
import sys
import csv
from datetime import datetime
csv.field_size_limit(sys.maxsize)


# def checkpam(argu, chr,start,secondaryhitpattern,strand, getgenome, querysequence, secondaryhitmismatch):
#     length = score = ggfound = 0
#     forwardseq = revcomplementseq = complementseq = "NONE"
#     for record in getgenome:
#         recordid = str(record.id)
#         #if (strand == 0 or strand <= 256) and recordid == chr:
#         if (strand == "+" ) and recordid == chr:
#             forwardseq = str(record.seq)[start-1:start+22]
#             if forwardseq.endswith(argu.pamsequence):
#                 ggfound = ggfound + 1
#             elif argu.pamsequence in forwardseq[-4:]:
#                 ggfound = ggfound + 1
#         elif (strand == "-") and recordid == chr:
#             complementseq = str(record.seq)[start-1:start + 22]
#             newseq = Seq(complementseq)
#             revcomplementseq = str(newseq.reverse_complement())
#             if revcomplementseq.endswith(argu.pamsequence):
#                 ggfound = ggfound + 1
#             elif argu.pamsequence in revcomplementseq[-4:]:
#                 ggfound = ggfound + 1
#     if ggfound > 0:
#         #print ("inside checkpam.. before getscore")
#         (length, score) = getscore(secondaryhitpattern, secondaryhitmismatch)
#     return(length, score, ggfound,forwardseq,revcomplementseq)

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

#
# def calculatewrite(argu, fileoutmap, gc_content, counthits, countoffhits, inwindscore, outwindscore, queryname, querychr, arbitraryid, querystart, offtargetlist, secondaryhitschr, bestwindow, newlist, querysequence, getstddev):
#     import string
#     suminwindscore = sum(map(float, inwindscore))
#     sumoutwindscore = sum(map(float, outwindscore))
#     finalscore = ((float(counthits) * (suminwindscore/float(counthits)))/float(counthits+countoffhits)) - ((float(countoffhits) * (sumoutwindscore/float(countoffhits)))/float(counthits+countoffhits))
#     ### Below sequence translate take more time ###
#     if counthits > argu.minhits:
#         if "candidaternars" in queryname:
#             ## For python 3.0 it changed to str.maketrans
#             #complements = string.maketrans('ATGC', 'TACG')
#             #complements = str.maketrans('ATGC', 'TACG')
#             #print (querysequence)
#             #querysequence = querysequence.translate(complements)[::-1]
#             #print (querysequence)
#             ############ Above sequence translate take more time ###
#             fileoutmap.write(str(queryname) + "_" + str(querychr) + "_" + str(arbitraryid) + "_" + str(querystart) + "\t" + str(querysequence) + "\t" + str(gc_content)+ "\t" + str(finalscore) + "\t" + str(getstddev) + "\t" + str(suminwindscore) + "\t" + str(sumoutwindscore) +  "\t" + str(counthits) + "\t" + str(len(offtargetlist)) + "\t" + str(querychr) + ":" + str(min(bestwindow)) + "-" + str(max(bestwindow)) + "\t" + str(newlist) + "\n")
#         elif "candidaternafs" in queryname:
#             fileoutmap.write(str(queryname) + "_" + str(querychr) + "_" + str(arbitraryid) + "_" + str(querystart) + "\t" + str(querysequence) + "\t" + str(gc_content) + "\t" + str(finalscore) + "\t" + str(getstddev) + "\t" + str(suminwindscore) + "\t" + str(sumoutwindscore) + "\t" + str(counthits) + "\t" + str(len(offtargetlist)) + "\t" + str(querychr) + ":" + str(min(bestwindow)) + "-" + str(max(bestwindow)) + "\t" + str(newlist) + "\n")
#     elif counthits < 1:
#         print ("Pattern Match error")

def calculatestdev(bestwindow):
    if len(bestwindow) > 4:
        getstddev = numpy.std(bestwindow, axis=0)
    else:
        getstddev = "NULL"
    return getstddev

def getallresults(argu):
    dataresult = pd.read_csv(os.path.join(argu.workingdirectory,"Crispr_broad_multihits.xls"), names=['Chromosome', 'End', 'Start', 'Totalhits', 'finalscore', 'offtargetshits', 'score', 'Sequence', 'stdev', 'stringforgroup', 'taggingcount'], sep="\t", header=None)
##Chromosome	End	Start	Totalhits	finalscore	offtargetshits	score	sequence	stdev	stringforgroup	taggingcount
    print ("Pooling the results...")
    column_names = ['Chromosome', 'Start', 'End', 'Totalhits', 'finalscore', 'offtargetshits', 'score', 'Sequence', 'stdev', 'stringforgroup', 'taggingcount']
    dataresult = dataresult.reindex(columns=column_names)
    print (dataresult)
    dataresult1 = dataresult.sort_values(['offtargetshits', 'stringforgroup'], ascending=[True, True])
    dataresult1.to_csv("Crispr_broad_multihits_sorted.xls", sep="\t", index=False, header=True)


def checkpandapamquerywindow (rowdatchunk, argu, splitinputfile):
    genomesplitoriginal = os.path.join(argu.workingdirectory,"Candidate_crisprnatabbed.txt")
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
        finalscore = ((float(counthits) * (suminwindscore / float(counthits))) / float(counthits + countoffhits)) - ((float(countoffhits) * (sumoutwindscore / float(countoffhits))) / float(counthits + countoffhits))
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
    genomesplitoriginal = os.path.join(argu.workingdirectory,"Candidate_crisprnatabbed.txt")
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
    ts = time.time()
    print(ts)
    print ("before overlap")
    overlappedggdataframe = (pr.PyRanges(eachrowdataframe3).overlap(pr.PyRanges(datachunkallcrna1))).df
    ts = time.time()
    print(ts)
    print ("After overlap")

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
        ##print (getscorewindow2)

        #getscorewindow1[['Chromosome', 'Start', 'End']] = getscorewindow1['stringforgroup'].str.split("-", expand=True)
        #print(getscorewindow2)
        ##########
        getscorenonwindow = (pr.PyRanges(overlappedggdataframe).join(filemakewindowsinwindow1, how="left")).df
        sumoutwindscore = getscorenonwindow[getscorenonwindow['Start_b'] < 0]['score'].sum()
        countoffhits = getscorenonwindow[getscorenonwindow['Start_b'] < 0].shape[0]
        #############

        #suminwindscore = sum(map(float, inwindscore))
        #sumoutwindscore = sum(map(float, outwindscore))
        finalscore = ((float(counthits) * (suminwindscore / float(counthits))) / float(counthits + countoffhits)) - ((float(countoffhits) * (sumoutwindscore / float(countoffhits))) / float(counthits + countoffhits))
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
    for datachunk in pd.read_csv(splitinputfile, header=None, names =['crnaid', 'strand', 'firsthitschr', 'firsthitstart', 'sequence', 'colCIGAR', 'Totalhits'], sep="\t", comment='@', dtype=str, engine='c', index_col=False, low_memory=False, na_filter=True, chunksize=200000, iterator=True):
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



#def checkpandapamsingle(rowdatchunk, argu, datachunkallcrna1):

def checkpandapamsingle(eachrowdataframe, argu, datachunkallcrna1):
#working module raj 25 april 2022
    #eachrowdataframe = pd.DataFrame(rowdatchunk['colCIGAR'], columns=['flips'])
    #eachrowdataframe['flips'] = eachrowdataframe['flips'].str.replace('XA:Z:', '', regex=False)
    eachrowdataframe['colCIGAR'] = eachrowdataframe['colCIGAR'].str.replace('+', '+,', regex=False)
    eachrowdataframe['colCIGAR'] = eachrowdataframe['colCIGAR'].str.replace('-', '-,', regex=False)
    eachrowdataframe['colCIGAR'] = eachrowdataframe['colCIGAR'].str.replace('XA:Z:', '', regex=False)
    #print(eachrowdataframe.head(10))
    print("Before split")
    eachrowdataframe[['Chromosome', 'Str', 'Start', 'mismatch', 'unusedstrings']] = eachrowdataframe['colCIGAR'].str.split(',', expand=True)
    if argu.avoidscoreoff == 1:
        eachrowdataframe['score'] = 1
        eachrowdataframe.loc[:, 'taggingcount'] = 1
    else:
        eachrowdataframe['score'] = eachrowdataframe.apply(lambda x: getscore(x, x['mismatch'], pd.to_numeric(x['unusedstrings']) ), axis=1)
        eachrowdataframe.loc[:,'taggingcount'] = 1
    #########################

    eachrowdataframe3 = eachrowdataframe[['Chromosome', 'Start', 'score','taggingcount','crnaid','sequence']]
    #print (eachrowdataframe3.head(10))
    eachrowdataframe3 = eachrowdataframe3.dropna()
    #eachrowdataframe3['End'] = eachrowdataframe3['Start'].astype(int) + 27
    eachrowdataframe3.loc[:,'End'] = eachrowdataframe3['Start'].astype(int) + 27
    #eachrowdataframe3 = eachrowdataframe3.dropna()
    #print (eachrowdataframe3.head(10))

    #datachunkallcrna1 = pd.DataFrame()
    #datachunkallcrna1['Chromosome'] = datachunkallcrna['Chromosome']
    #datachunkallcrna1['Start'] = datachunkallcrna['Startggposition']
    #datachunkallcrna1.loc[:,'Start'] = datachunkallcrna1['Start'].astype(int) - 4
    #datachunkallcrna1.loc[:,'End'] = datachunkallcrna['Startggposition'].astype(int) + 4
    #now = datetime.now()
    #current_time = now.strftime("%H:%M:%S")
    #print(current_time)
    #print ("before overlap")
    overlappedggdataframe = (pr.PyRanges(eachrowdataframe3).overlap(pr.PyRanges(datachunkallcrna1))).df
    #print ("overlap done")
    #now = datetime.now()
    #current_time = now.strftime("%H:%M:%S")
    #print(current_time)
    #print(overlappedggdataframe.head(10))
    foreachoverlappedggdataframe = overlappedggdataframe[['Chromosome','Start', 'End', 'score','taggingcount','crnaid','sequence']].groupby('crnaid')
    #print(foreachoverlappedggdataframe.head(10))
    for eachname, eachgroupedframe in foreachoverlappedggdataframe:
        #print(eachname)
        #print (eachgroupedframe.head(10))
        if (len(eachgroupedframe) >= int(argu.minhits)):
            allchrmax = eachgroupedframe[['Chromosome', 'End']].groupby('Chromosome').End.agg('max').reset_index()
            allchrmin = eachgroupedframe[['Chromosome', 'Start']].groupby('Chromosome').Start.agg('min')
            allchrnew = pd.merge(allchrmin, allchrmax, on='Chromosome')
            allchrnew.columns = ['Chromosome', 'Start', 'End']
            allchrnew = allchrnew[~allchrnew.isin([np.nan, np.inf, -np.inf]).any(1)]
            filemakewindowsinwindow1 = pr.PyRanges(allchrnew).window(argu.windowsize)
            if argu.avoidscoreoff == 1:
                getscorewindow = (filemakewindowsinwindow1.join(pr.PyRanges(eachgroupedframe))).df
                #print(getscorewindow.head(10))
                getscorewindow['stringforgroup'] = getscorewindow[['Chromosome', 'Start', 'End']].astype(str).agg('-'.join, axis=1)
                getscorewindowsum = getscorewindow.groupby(['stringforgroup'], as_index=False)['score'].sum()
                getscorewindowcount = getscorewindow.groupby(['stringforgroup'], as_index=False)['taggingcount'].count()
                getscorewindow2 = pd.merge(getscorewindowcount, getscorewindowsum, on='stringforgroup')
                #print(getscorewindow2.head(10))
                #getscorewindow['stringforgroup'] = getscorewindow[['Chromosome', 'Start', 'End']].astype(str).agg('-'.join, axis=1)
                #getscorewindow2.loc[:,"sequence1"] = eachgroupedframe['sequence'][0]
                getscorewindow2.loc[:,"Totalhits"] = len(eachgroupedframe)
                getscorewindow2.loc[:,"sequence"] = eachgroupedframe['sequence'].iloc[0]
                #print(hitslistforstdevdf)
                hitslistforstdevdf = getscorewindow.groupby('stringforgroup')['Start_b'].apply(list).reset_index(name='starthistaslist')
                getscorewindow2["stdev"] = hitslistforstdevdf.apply(lambda x: calculatestdev( x['starthistaslist']), axis=1)
                getscorewindow2.loc[:,"finalscore"] = 1 
                #argu.windownumbers = 1
                #print(getscorewindow2.head(10))
                #print(eachgroupedframe['sequence'].iloc[0])
                #getscorewindow['stringforgroup'] = getscorewindow[['Chromosome', 'Start', 'End']].astype(str).agg('-'.join, axis=1)
                getscorewindow3 = getscorewindow2.sort_values(['taggingcount'], ascending=[False]).head(argu.windownumbers)
                getscorewindow3['offtargetshits'] = pd.to_numeric(getscorewindow3['Totalhits']).iloc[0] - getscorewindow3['taggingcount'].sum()
                getscorewindow3[['Chromosome', 'Start', 'End']] = getscorewindow3['stringforgroup'].str.split("-", expand=True)
                fileoutstd = open(os.path.join(argu.workingdirectory, "Crispr_broad_multihits.xls"), 'a')
                getscorewindow3 = getscorewindow3.reindex(sorted(getscorewindow3.columns), axis=1)
                getscorewindow3.to_csv(fileoutstd, sep="\t", index=False, header=False)
                #print (getscorewindow3.head(10))
                #now = datetime.now()
                #current_time = now.strftime("%H:%M:%S")
                #print(current_time)
            else:
                getscorewindow = (filemakewindowsinwindow1.join(pr.PyRanges(eachgroupedframe))).df
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
                ########## OFF Target
                getscorenonwindow = (pr.PyRanges(eachgroupedframe).join(filemakewindowsinwindow1, how="left")).df
                sumoutwindscore = getscorenonwindow[getscorenonwindow['Start_b'] < 0]['score'].sum()
                countoffhits = getscorenonwindow[getscorenonwindow['Start_b'] < 0].shape[0]
                #############

                #suminwindscore = sum(map(float, inwindscore))
                #sumoutwindscore = sum(map(float, outwindscore))
                finalscore = ((float(counthits) * (suminwindscore / float(counthits))) / float(counthits + countoffhits)) - ( (float(countoffhits) * (sumoutwindscore / float(countoffhits))) / float(counthits + countoffhits))
                ###########
                #print (finalscore)
                #getstddev = calculatestdev(bestwindow)
                getscorewindow2.loc[:,"finalscore"] =  finalscore

                #getscorewindow2["crnaid"] = eachgroupedframe['crnaid']
                getscorewindow2.loc[:,"Totalhits"] = len(eachgroupedframe)
                getscorewindow2.loc[:,"sequence"] = eachgroupedframe['sequence'].iloc[0]
                #argu.windownumbers = 1
                #print(getscorewindow2)
                ####################
                argu.windownumbers = 1
                getscorewindow3 = getscorewindow2.sort_values(['taggingcount'], ascending=[False]).head(argu.windownumbers)
                # print(data1['Totalhits'].iloc[0])
                # print(data1)
                getscorewindow3['offtargetshits'] = pd.to_numeric(getscorewindow3['Totalhits']).iloc[0] - getscorewindow3['taggingcount'].sum()
                # data1['off-targets-hits'] = data1.apply(lambda row: TotalHitsforthatRNA - inwindowtarget, axis=1)
                # print ("I am writing")
                getscorewindow3[['Chromosome', 'Start', 'End']] = getscorewindow3['stringforgroup'].str.split("-", expand=True)
                #print (getscorewindow3)
                fileoutstd = open(os.path.join(argu.workingdirectory, "Crispr_broad_multihits.xls"), 'a')
                getscorewindow3 = getscorewindow3.reindex(sorted(getscorewindow3.columns), axis=1)
                getscorewindow3.to_csv(fileoutstd, sep="\t", index=False, header=False)
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                print(current_time)





def overlapeachchromosomesingle(argu, splitinputfile):
    print ("Processing:",splitinputfile)
    genomesplitoriginal = os.path.join(argu.workingdirectory,"Candidate_crisprnatabbed.txt2col")
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(current_time)
    ###datachunkallcrna = pd.read_csv(genomesplitoriginal, header=None, names=['name', 'Str', 'Chromosome', 'crnanumber', 'Startggposition'], sep="\t", dtype=str, engine='c', index_col=False, low_memory=False, na_filter=True)
    datachunkallcrna = pd.read_csv(genomesplitoriginal, header=None, names=['Chromosome', 'Startggposition'], sep="\t", dtype=str, engine='c', index_col=False, low_memory=False, na_filter=True)
    print ("Two col data loaded")
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(current_time)
    datachunkallcrna1 = pd.DataFrame()
    datachunkallcrna1['Chromosome'] = datachunkallcrna['Chromosome']
    datachunkallcrna1['Start'] = datachunkallcrna['Startggposition']
    datachunkallcrna1.loc[:,'Start'] = datachunkallcrna1['Start'].astype(int) - 4
    datachunkallcrna1.loc[:,'End'] = datachunkallcrna['Startggposition'].astype(int) + 4
    for datachunk in pd.read_csv(splitinputfile, header=None, usecols =['crnaid','colCIGAR','sequence'], names =['crnaid', 'strand', 'firsthitschr', 'firsthitstart', 'sequence','waste1','waste2', 'colCIGAR', 'Totalhits'], sep="\t", comment='@', dtype=str, engine='c', index_col=False, low_memory=False, na_filter=True, chunksize=2000, iterator=True):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print(current_time)
        print ("Filtered hits started after loading chunk1")
        #datachunk["colCIGAR"] = datachunk["colCIGAR"].str.split(";", expand=False)
        datachunk['colCIGAR'] = datachunk["colCIGAR"].str.split(";")
        datachunk = datachunk[datachunk["colCIGAR"].str.len() > int(argu.minhits)]
        datachunk = datachunk.set_index(['crnaid','sequence']).apply(pd.Series.explode).reset_index()
        datachunk1 = datachunk[['crnaid', 'colCIGAR','sequence']]
        #print (datachunk1.head(10))
        #datachunk = datachunk[datachunk["colCIGAR"].str.len() > int(argu.minhits)]
        #checkpandapamsingle(datachunk1, argu, datachunkallcrna1)
        checkpandapamsingle(datachunk1, argu, datachunkallcrna1)
        ###datachunk.apply(lambda x: checkpandapamsingle(x, argu, splitinputfile, datachunkallcrna1), axis =1 )



def getmultisgrna(argu):
    #multisgrnafile = "crispr_broad_multisgrna.xls"
    multisgrnafile= pr.PyRanges(pd.read_csv(argu.crisprbroadresfile, sep="\t"))
    getscorewindow = (multisgrnafile.join(multisgrnafile,slack=-10)).df
    #getscorewindow1 = getscorewindow.groupby(['crnaid', 'crnaid_b']).agg(['count'])
    getscorewindow1 = getscorewindow.groupby(['crnaid', 'crnaid_b'])['crnaid'].count().reset_index(name='counts').sort_values(['counts'], ascending=False)
    #getscorewindow2 = getscorewindow1.sort_values(['count'], ascending=[False]).head(argu.windownumbers)
    #getscorewindow1 = getscorewindow.groupby(['crnaid', 'crnaid_b'])['crnaid'].count().reset_index(name='counts')
    getscorewindow2 = getscorewindow1.groupby(['crnaid']).head(argu.sgrnanumbers)
    getscorewindow3 = pd.merge(getscorewindow2, getscorewindow, how='left', left_on=['crnaid','crnaid_b'], right_on = ['crnaid','crnaid_b'] )
    #overlappedggdataframe = (pr.PyRanges(eachrowdataframe2).overlap(pr.PyRanges(datachunkallcrna1))).df
    #getscorewindow.to_csv(os.path.join(argu.workingdirectory, "crispr_broad_multisgrna.xls"), sep="\t", index=False, header=True)
    #getscorewindow1.to_csv(os.path.join(argu.workingdirectory, "crispr_broad_grouped_multisgrna.xls"), sep="\t", index=False,header=True)
    #getscorewindow2.to_csv(os.path.join(argu.workingdirectory, "crispr_broad_grouped_multisgrna2.xls"), sep="\t", index=False, header=True)
    getscorewindow3.to_csv(os.path.join(argu.workingdirectory, "Crispr_broad_multisgrna.xls"), sep="\t", index=False, header=True)
    #filemultisgrnain = open(os.path.join(argu.workingdirectory, "crispr_broad_multisgrna.xls"), 'w+')
    #p2 = subprocess.Popen(['bedtools', 'intersect', '-a', str(argu.crisprbroadresfile), '-b', str(argu.crisprbroadresfile), '-wao'], stdout=filemultisgrnain)
    #p2.communicate()







def checkpandapamsingleWORKINGOLD(rowdatchunk, argu, splitinputfile, datachunkallcrna1):
    #genomesplitoriginal = "Candidate_crisprnatabbed.txt"
    #datachunkallcrna = pd.read_csv(genomesplitoriginal, header=None, names=['name', 'Str', 'Chromosome', 'crnanumber', 'Startggposition'], sep="\t", dtype=str, engine='c', index_col=False, low_memory=False, na_filter=True)
    eachrowdataframe = pd.DataFrame(rowdatchunk['colCIGAR'], columns=['flips'])
    eachrowdataframe['flips'] = eachrowdataframe['flips'].str.replace('XA:Z:', '', regex=False)
    eachrowdataframe['flips'] = eachrowdataframe['flips'].str.replace('+', '+,', regex=False)
    eachrowdataframe['flips'] = eachrowdataframe['flips'].str.replace('-', '-,', regex=False)
    eachrowdataframe = eachrowdataframe['flips'].str.replace('XA:Z:', '', regex=False)
    eachrowdataframe1 = pd.DataFrame(eachrowdataframe.str.split(',').tolist(), columns=['Chromosome', 'Str', 'Start', 'mismatch', 'unusedstrings'])

    eachrowdataframe2 = eachrowdataframe1.copy()
    eachrowdataframe2.loc[eachrowdataframe1.index.max()+1] = [rowdatchunk['chr'],rowdatchunk['strand'], rowdatchunk['firsthitstart'], '23M','0']
    eachrowdataframe2 = eachrowdataframe2.dropna()
    if argu.avoidscoreoff == 1:
        eachrowdataframe2['score'] = 1
        eachrowdataframe2.loc[:, 'taggingcount'] = 1
    else:
        eachrowdataframe2['score'] = eachrowdataframe2.apply(lambda x: getscore(x, x['mismatch'], pd.to_numeric(x['unusedstrings']) ), axis=1)
        eachrowdataframe2.loc[:,'taggingcount'] = 1
    #########################

    eachrowdataframe3 = eachrowdataframe2[['Chromosome', 'Start', 'score','taggingcount']]
    #eachrowdataframe3['End'] = eachrowdataframe3['Start'].astype(int) + 27
    eachrowdataframe3.loc[:,'End'] = eachrowdataframe3['Start'].astype(int) + 27
    eachrowdataframe3 = eachrowdataframe3.dropna()
    #print (eachrowdataframe3)

    #datachunkallcrna1 = pd.DataFrame()
    #datachunkallcrna1['Chromosome'] = datachunkallcrna['Chromosome']
    #datachunkallcrna1['Start'] = datachunkallcrna['Startggposition']
    #datachunkallcrna1.loc[:,'Start'] = datachunkallcrna1['Start'].astype(int) - 4
    #datachunkallcrna1.loc[:,'End'] = datachunkallcrna['Startggposition'].astype(int) + 4
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(current_time)
    print ("before overlap")
    overlappedggdataframe = (pr.PyRanges(eachrowdataframe3).overlap(pr.PyRanges(datachunkallcrna1))).df
    print ("overlap done")
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(current_time)
    #print (overlappedggdataframe)
    if (len(overlappedggdataframe) >= int(argu.minhits)):
        allchrmax = overlappedggdataframe[['Chromosome', 'End']].groupby('Chromosome').End.agg('max').reset_index()
        allchrmin = overlappedggdataframe[['Chromosome', 'Start']].groupby('Chromosome').Start.agg('min')

        allchrnew = pd.merge(allchrmin, allchrmax, on='Chromosome')
        allchrnew.columns = ['Chromosome', 'Start', 'End']
        allchrnew = allchrnew[~allchrnew.isin([np.nan, np.inf, -np.inf]).any(1)]
        filemakewindowsinwindow1 = pr.PyRanges(allchrnew).window(argu.windowsize)
        #print (filemakewindowsinwindow1)
        #filemakewindowsinwindow1 = pr.PyRanges(filemakewindowsinwindow)

        getscorewindow = (filemakewindowsinwindow1.join(pr.PyRanges(overlappedggdataframe))).df
        #print (getscorewindow)
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
        finalscore = ((float(counthits) * (suminwindscore / float(counthits))) / float(counthits + countoffhits)) - ( (float(countoffhits) * (sumoutwindscore / float(countoffhits))) / float(counthits + countoffhits))
        ###########
        #print (finalscore)
        #getstddev = calculatestdev(bestwindow)
        getscorewindow2.loc[:,"finalscore"] =  finalscore

        getscorewindow2["crnaid"] = rowdatchunk['crnaid']
        getscorewindow2["sequence"] = rowdatchunk['sequence']
        getscorewindow2.loc[:,"Totalhits"] = len(overlappedggdataframe)
        #print(getscorewindow2)
        ####################
        argu.windownumbers = 1
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


