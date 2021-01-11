
def readdirectories(subparsers):
    """ Read bed/bam files from a directory. """
    """ Options saved in argu. """

    parser_t = subparsers.add_parser('genomesplit', help="use python epicrisprtarget.py genomesplit -h")
    parser_t.add_argument('-f', '--genome_fasta', required=True, type=str, dest='genomesplitfasta', action='store', help="Genome sequence in FASTA format")
    parser_t.add_argument('-c', '--chromosome', required=False, type=str, dest='chr', action='store', help="chromosome name")
    parser_t.add_argument('-d', '--working_directory', required=True, type=str, dest='workingdirectory', action='store',help="Complete path of output directory")
    parser_t.add_argument('-p', '--pam_sequence', required=True, type=str, dest='pamsequence', action='store', help="PAM sequence string eg: GG for NGG PAM")
    parser_t.add_argument('-l', '--get_candidaternalength', required=True, type=int, dest='candidaternalength', action='store', help="Candidate gRNA length")
    parser_t.add_argument('-g', '--get_gc', required=True, type=int, default=45, dest='gc', action='store',help="gc content in integer")
    parser_t.add_argument('-t', '--num_threads', required=False, type=int, default=1, dest='threads', action='store', help="Launch t number of threads in parallel")


    parser_c = subparsers.add_parser('createindex', help="use python epicrisprtarget.py createindex -h")
    parser_c.add_argument('-f', '--genome_fasta', required=True, type=str, dest='genomesplitfasta', action='store', help="Genome sequence in FASTA format")

    parser_f = subparsers.add_parser('maptogenome', help="use python epicrisprtarget.py maptogenome -h")
    parser_f.add_argument('-f', '--genome_fasta', required=True, type=str, dest='genomesplitfasta', action='store', help="Genome sequence in FASTA format")
    parser_f.add_argument('-d', '--working_directory', required=True, type=str, dest='workingdirectory', action='store',help="Complete path of output directory")
    parser_f.add_argument('-m', '--num_mismatch', required=True, type=int, dest='mismatch', action='store',help="Maximum number of mismatches for Bowtie2 alignment")
    parser_f.add_argument('-t', '--num_threads', required=True, type=int, dest='threads', action='store',help="Launch t number of threads in parallel")
    #parser_f.add_argument('-s', '--get_sequence', required=True, type=str, dest='getsequences', action='store',help="Get the mismatched hit sequence")
    parser_f.add_argument('-n', '--get_minhits', required=True, type=int, dest='minhits', action='store',help="Minimum number of hits in a window (filter)")
    parser_f.add_argument('-k', '--get_maxmum_total_hits', required=True, type=int, default=10000, dest='maxhits', action='store', help="Maximum total number of hits")
    parser_f.add_argument('-g', '--get_gc', required=True, type=int, dest='gc', action='store',help="gc content in integer")
    parser_f.add_argument('-l', '--get_candidaternalength', required=True, type=int, dest='candidaternalength', action='store', help="Candidate gRNA length")

    parser_f = subparsers.add_parser('filterhits', help="use python epicrisprtarget.py filterhits -h")
    parser_f.add_argument('-d', '--working_directory', required=True, type=str, dest='workingdirectory', action='store',help="Complete path of output directory")
    parser_f.add_argument('-t', '--num_threads', required=True, type=int, dest='threads', action='store',help="Launch t number of threads in parallel")
    parser_f.add_argument('-n', '--get_minhits', required=True, type=int, dest='minhits', action='store',help="Minimum number of hits in a window (filter)")

    parser_t = subparsers.add_parser('findwindow', help="use python epicrisprtarget.py findwindow -h")
    parser_t.add_argument('-f', '--genome_fasta', required=True, type=str, dest='genomesplitfasta', action='store', help="Genome sequence in FASTA format")
    parser_t.add_argument('-d', '--working_directory', required=True, type=str, dest='workingdirectory', action='store',help="Complete path of output directory")
    parser_t.add_argument('-p', '--pam_sequence', required=True, type=str, dest='pamsequence', action='store', help="PAM sequence string eg: NGG")
    parser_t.add_argument('-t', '--num_threads', required=True, type=int, dest='threads', action='store',help="Launch t number of threads in parallel")
    parser_t.add_argument('-n', '--get_minhits', required=True, type=int, dest='minhits', action='store',help="Minimum number of hits in a window (filter)")
    parser_t.add_argument('-w', '--get_window', required=True, type=int, dest='window', action='store',help="Window size in bp")
    parser_t.add_argument('-l', '--get_candidaternalength', required=True, type=int, dest='candidaternalength', action='store', help="Candidate gRNA length")

def arg_parsing():
    """ Parse the subcommand along with its arguments. """
    descr = '''
    Finds gRNA sequences in large windows.
    '''
    import argparse
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers( title='Subcommands', dest="subcommand")
    readdirectories(subparsers)
    # second_parsing(subparsers)
    # third_parsing(subparsers)
    argu = parser.parse_args()
    return argu