#! /usr/bin/python
from argumentparse import *
from mapbwa import *
from readfasta import *
from checkpamandscore import *
import threading
import time
import multiprocessing
from multiprocessing import Process, Queue
from functools import partial

##############################################################################
#                               MAIN
##############################################################################
if __name__ == "__main__":
    start = time.time()
    argu = arg_parsing()
    # ####### Check all here #####
    # ####### Multithreading here #####
    if argu.subcommand == 'genomesplit':
        if argu.inputbed is not None:
            #print (argu.readbedtofasta#)
            convertquerybedtofasta = readbedtofasta(argu)
            counttotalgg = createtabbedpamfile(argu)
            countinputlines = readgenomefastafile(argu, convertquerybedtofasta)
            splitinputfiles = splitinputfileformulti(argu,countinputlines)
        if argu.inputbed is None:
            countinputlines = readgenomefastafile(argu, argu.genomesplitfasta)
            counttotalgg = createtabbedpamfile(argu)
            splitinputfiles = splitinputfileformulti(argu,countinputlines)
    if argu.subcommand == 'createindex':
        indexfastatogenome(argu)
    if argu.subcommand == 'maptogenome':
        splitinputfiles = getlistfiles(argu)
        print (splitinputfiles)
        newthreads = []
        for splitinputfile in splitinputfiles:
            newthread = Process(target=mapbwafastatogenome, args=(argu, splitinputfile))
            newthreads.append(newthread)
        for newthread in newthreads:
            newthread.start()
        for newthread in newthreads:
            newthread.join()
    if argu.subcommand == 'filterhits':
        getminhitspattern = getminhits(argu)
        splitinputfiles = getsamlistfiles(argu)
        # pool = multiprocessing.Pool(processes=40)
        # newpoolfunc = partial(filtersam, argu, getminhitspattern)
        # newpool = pool.map(newpoolfunc, splitinputfiles)
        #print (splitinputfiles)
        newthreads = []
        for splitinputfile in splitinputfiles:
            newthread = Process(target=filtersam, args=(argu, splitinputfile))
            newthread.Daemon = True
            newthreads.append(newthread)
        for newthread in newthreads:
            newthread.start()
        for newthread in newthreads:
            newthread.join()
        #print (splitinputfiles)
        #newthreads = []
        # for splitinputfile in splitinputfiles:
        #     newthread = threading.Thread(target=filtersam, args=(argu, splitinputfile, getminhitspattern))
        #     newthreads.append(newthread)

        #newpool = pool.map(overlapeachchromosome, splitinputfiles)
        # for newthread in newthreads:
        #     newthread.start()
        # for newthread in newthreads:
        #     newthread.join()
    if argu.subcommand == 'findwindow':
        newthreads = []
        splitinputfiles = getsamfilteredfiles(argu)
        for splitinputfile in splitinputfiles:
            newthread = Process(target=overlapeachchromosomesingle, args=(argu, splitinputfile))
            newthreads.append(newthread)
        for newthread in newthreads:
            newthread.daemon = True
            newthread.start()
        for newthread in newthreads:
            newthread.join()
        getallresults(argu)

    if argu.subcommand == 'findmultiwindow':
        newthreads = []
        splitinputfiles = getsamfilteredfiles(argu)
        for splitinputfile in splitinputfiles:
            newthread = Process(target=overlapeachchrommultiwindow, args=(argu, splitinputfile))
            newthreads.append(newthread)
        for newthread in newthreads:
            newthread.daemon = True
            newthread.start()
        for newthread in newthreads:
            newthread.join()
        getallresults(argu)
    if argu.subcommand == 'multisgrna':
        getmultisgrna(argu)
        #newthread1.start()
        #newthread1.join()
        # newthreads = []
        # for splitinputfile in splitinputfiles:
        #     newthread = Process(target=multisgrna, args=(argu))
        #     newthreads.append(newthread)
        # for newthread in newthreads:
        #     newthread.daemon = True
        #     newthread.start()
        # for newthread in newthreads:
        #     newthread.join()
        # getallresults(argu)
    end = time.time()
    print ("Job Finished:", end-start, "seconds")

