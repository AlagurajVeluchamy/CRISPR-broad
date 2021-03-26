from subprocess import call
import pandas as pd
import sys
import csv
csv.field_size_limit(sys.maxsize)

def indexfastatogenome(argu):
    ##Indexing the genome
    print ("BWA indexing the genome..")
    #call(['bowtie2-build', argu.genomesplitfasta, argu.genomesplitfasta])
    call(['bwa', 'index', argu.genomesplitfasta])

def mapbwafastatogenome(argu,splitinputfile):
    splitinputfilesam = splitinputfile + ".sam"
    splitinputfilealn = splitinputfile + ".aln"
    ##Mapping cgRNA to the genome
    print ("BWA based mapping of ",splitinputfile," to the genome..")
    call(['bwa', 'aln', '-t', str(argu.threads), '-n', str(argu.mismatch), '-k', str(argu.mismatch), '-l', str(argu.candidaternalength), '-N', '-f', splitinputfilealn, argu.genomesplitfasta, splitinputfile])
    call(['bwa', 'samse', '-n', str(argu.maxhits), '-f', splitinputfilesam, argu.genomesplitfasta, splitinputfilealn, splitinputfile])
    ### bowtiecall(['bowtie2', '-x', argu.genomesplitfasta, '-N', '1','-f', splitinputfile, '--xeq', '-a', '-S', str(splitinputfile+".sam")])


def filtersambowtie(argu, splitinputfile,getminhitspattern):
    colnames=['crnaid','strand', 'chr', 'hitstart', 'waste2','cigar','waste3', 'waste4', 'waste5', 'sequence']
    data = pd.read_csv(splitinputfile, sep="\t", names=colnames,usecols=['crnaid', 'strand', 'chr', 'hitstart', 'waste2', 'cigar', 'waste3', 'waste4', 'waste5', 'sequence'], skiprows = 8, comment='@', dtype = str)
    data['merged'] = data['chr'].astype(str) + '_' + data['hitstart'].astype(str) + '_' + data['cigar'].astype(str) + '_' + data['strand'].astype(str)
    data1 = data[['crnaid', 'merged', 'sequence']].groupby('crnaid').agg({'sequence': 'first', 'merged': lambda x: ','.join(x)}).reset_index()[['crnaid', 'sequence', 'merged']]
    datafiltered = data1[data1['merged'].str.contains(getminhitspattern, na=False)]
    datafiltered.to_csv(splitinputfile+"_filtered.txt", sep="\t", index=False, header=False)

def filtersam(argu, splitinputfile):
    print ("filtering for number of hits:", splitinputfile)
    colnames = ['crnaid', 'strand', 'chr', 'hitstart', 'col5', 'col6', 'col7', 'col8', 'col9', 'sequence',
                'col11', 'col12', 'col13', 'col14', 'col15', 'col16', 'col17', 'col18',
                'col19', 'allhits']

    #data = pd.read_csv(splitinputfile, header=None, usecols=['crnaid', 'strand', 'chr', 'allhits', 'sequence', 'col14'],
    #    names=colnames, sep="\t", skiprows=8, comment='@', dtype=str, engine='c', index_col=False,low_memory=False, na_filter=False)
    # datasplit = data["col14"].str.split("X0:i:", n = 1, expand = True)
    # data["counthits"] = datasplit[1].apply(pd.to_numeric)
    # datafiltered = data.loc[data['counthits'] >= 5]
    # print (datafiltered.head(5))
    # datafiltered.to_csv(splitinputfile + "_filtered.txt", sep="\t", index=False, header=False)
    for datachunk in pd.read_csv(splitinputfile, header=None, usecols=['crnaid', 'strand', 'chr', 'hitstart', 'allhits', 'sequence', 'col14', 'col15'],
        names=colnames, sep="\t", comment='@', dtype=str, engine='c', index_col=False, low_memory=False,
        na_filter=False, chunksize=200000, iterator=True):
        datasplit1 = datachunk["col14"].str.split("X0:i:", n = 1, expand = True)
        if datasplit1 is None:
            print(splitinputfile, "notworking")
        datasplit2 = datachunk["col15"].str.split("X1:i:", n=1, expand=True)
        datachunk["counthitopt"] = datasplit1[1].apply(pd.to_numeric)
        datachunk["counthitsubopt"] = datasplit2[1].apply(pd.to_numeric)
        #print (argu.minhits)
        #datafiltered = datachunk.loc[(datachunk['counthitopt'] + datachunk['counthitsubopt'] + 1)  >= int(argu.minhits)]
        datafiltered = datachunk[(((datachunk['counthitopt'] + datachunk['counthitsubopt'] + 1) >= int(argu.minhits)) & ((datachunk['counthitopt'] + datachunk['counthitsubopt'] + 1) <= int(argu.maxhits)))]
        datafiltered.to_csv(splitinputfile + "_filtered.txt", sep="\t", index=False, header=False, mode='a')
