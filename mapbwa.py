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
                'col19', 'CIGAR']

    #data = pd.read_csv(splitinputfile, header=None, usecols=['crnaid', 'strand', 'chr', 'allhits', 'sequence', 'col14'],
    #    names=colnames, sep="\t", skiprows=8, comment='@', dtype=str, engine='c', index_col=False,low_memory=False, na_filter=False)
    # datasplit = data["col14"].str.split("X0:i:", n = 1, expand = True)
    # data["counthits"] = datasplit[1].apply(pd.to_numeric)
    # datafiltered = data.loc[data['counthits'] >= 5]
    # print (datafiltered.head(5))
    # datafiltered.to_csv(splitinputfile + "_filtered.txt", sep="\t", index=False, header=False)
    for datachunk in pd.read_csv(splitinputfile, header=None, usecols=['crnaid', 'strand', 'chr', 'hitstart', 'CIGAR', 'sequence', 'col14', 'col15'], names=colnames, sep="\t", comment='@', dtype=str, engine='c', index_col=False, low_memory=False, na_filter=False, chunksize=200000, iterator=True):
        #datasplit1 = datachunk["col14"].str.split("X0:i:", n = 1, expand = True)
        #print (datachunk.head(10))
        #if datasplit1 is None:
        #print(splitinputfile, "notworking")
        #datasplit2 = datachunk["col15"].str.split("X1:i:", n=1, expand=True)
        #datachunk["counthitopt"] = datasplit1[1].apply(pd.to_numeric)
        #datachunk["counthitsubopt"] = datasplit2[1].apply(pd.to_numeric)
        datachunk['crnaidchronly'] = datachunk['crnaid'].str.split('_').str[1]
        #print(datachunk[["allhits","crnaidchronly"]])
        #datachunk['samechrhitscount'] = datachunk.text.str.count('happy')
        #datachunk['samechrhitscount'] = datachunk[["col15","crnaidchronly"]].apply(lambda x: x[0].count(f";{x[1]}"), axis=1)
        datachunk['samechrhitscount'] = datachunk[["CIGAR","crnaidchronly"]].apply(lambda x: x[0].count(f";{x[1]},"), axis=1)
        datachunk['samechrhitscount1'] = datachunk[["CIGAR","crnaidchronly"]].apply(lambda x: x[0].count(f":{x[1]},"), axis=1)
        datachunk['totalmultihitscount'] = datachunk[["CIGAR","crnaidchronly"]].apply(lambda x: x[0].count(f";"), axis=1) + 1
        datahitsfil = datachunk[datachunk["crnaidchronly"]  == datachunk["chr"]]
        datahitsfilt = datahitsfil[(datahitsfil["samechrhitscount"] + datahitsfil["samechrhitscount1"] + 1) == (datahitsfil["totalmultihitscount"])]
        #print (datachunk)
        #datafiltered = datachunk.loc[(datachunk['counthitopt'] + datachunk['counthitsubopt'] + 1)  >= int(argu.minhits)]
        datafiltered = datahitsfilt[((datahitsfilt['totalmultihitscount']  >= int(argu.minhits)) & (datahitsfilt['totalmultihitscount'] <= int(argu.maxhits)))]
        datafiltered1 = datafiltered[['crnaid', 'strand', 'chr', 'hitstart', 'sequence', 'col14', 'col15','CIGAR', 'totalmultihitscount']]
        datafiltered1.to_csv(splitinputfile + "_filtered.txt", sep="\t", index=False, header=False, mode='a')
