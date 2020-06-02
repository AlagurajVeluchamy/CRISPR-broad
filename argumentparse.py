from readfasta import *

def readdirectories(subparsers):
    """ Read bed/bam files from a directory. """
    """ Options saved in argu. """
    help_str = "Read bed/bam files from a directory."
    parser_t = subparsers.add_parser('genomesplit', help=help_str)
    parser_t.add_argument('-d', '--genome_fasta', required=False, type=str, dest='genomesplitfasta', action='store', help=help_str)
    parser_t.add_argument('-p', '--pam_sequence', required=False, type=str, dest='pamsequence', action='store', help=help_str)
    parser_t.add_argument('-o', '--output_directory', required=False, type=str, dest='outputdir', action='store',help=help_str)
    parser_t.add_argument('-m', '--num_mismatch', required=False, type=str, dest='mismatch', action='store',help=help_str)
    parser_t.add_argument('-t', '--num_threads', required=False, type=str, dest='threads', action='store',help=help_str)
    parser_t.add_argument('-s', '--get_sequence', required=False, type=str, dest='getsequences', action='store',help=help_str)
    parser_t.add_argument('-n', '--get_minhits', required=False, type=int, dest='minhits', action='store',help=help_str)
    parser_t.add_argument('-w', '--get_window', required=False, type=int, dest='window', action='store',help=help_str)
    parser_t.add_argument('-l', '--get_candidaternalength', required=False, type=int, dest='candidaternalength', action='store', help=help_str)


    ##parser_t.add_argument('-b', '--directory_bamfiles', required=False, type=str, dest='dir_bam', action='store', help=help_str)
    parser_t.set_defaults(func=readgenomefastafile)
def arg_parsing():
    """ Parse the subcommand along with its arguments. """
    descr = '''
    Find complexes from ChIP-seq mapped files.
    '''
    import argparse
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers( title='Subcommands')
    #parser.add_argument("square", type=int, help="display the square of a given number")
    readdirectories(subparsers)
    # tffm_apply_arg_parsing(subparsers)
    # pssm_train_arg_parsing(subparsers)
    # pssm_apply_arg_parsing(subparsers)
    # binary_train_arg_parsing(subparsers)
    # binary_apply_arg_parsing(subparsers)
    argu = parser.parse_args()
    #print argu
    return argu