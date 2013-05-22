"""Wrapper for BFAST.

Kenny Daily, 2012

"""

import subprocess
import argparse
import os
import Task

from ccrngspy import utils

logger = utils.make_local_logger("BFAST logging", level="debug", color="green")

def setup_scratch(fxn):
    """Wrapper to create temp directory.

    """

    def wrapper(self, args):
        try:
            os.mkdir(args.bfast_temp_dir)
        except OSError:
            pass
            
        return fxn(self, args)

    return wrapper
                
class BFASTBase(Task.Task):
    """Container for BFAST tasks.
    
    """
    
    _cmd = "/usr/local/bfast-0.7.0a/bin/bfast"
    
    def __init__(self, args=None):
        if args:
            self.args = args
        # reads_file=None, output_file=None, threads=1, space=0, main_indexes=None, other_params=None,
        # match_file=None, align_file_name=None, align_algorithm=0, pairing_strandedness=0, pair_positioning=0,
        # pairing_options=None, scoring_matrix_file=None, insert_size_avg=None, insert_size_stdev=None):


    def set_exec(self, cmd):
        """Change the executable to use.

        """

        self._cmd = cmd

    def argparse(self, parser):
        """Add BFAST match option group to an OptionParser.
        
        """
        parser.add_argument("--bfast_reference_fasta", type=str, default=None,
                            help="FASTA reference file (needs prebuilt index; -f from bfast).")
        
        parser.add_argument("--bfast_threads", type=int, default=1,
                            help="Specifies the number of threads (-n from bfast).")

        parser.add_argument("--bfast_output_file", type=str, help="BFAST output file.", default=None)

        parser.add_argument("--bfast_temp_dir", type=str, help="BFAST temporary file directory.", default=None)

        parser.add_argument("--bfast_other_params", dest="bfast_other_params", type=str, default=None,
                            help="Other bfast params.")

        subparsers = parser.add_subparsers(help="BFAST sub-command help", dest='subparser_name')

        matchparser = subparsers.add_parser("match", help="match help")

        matchparser.add_argument("--bfast_reads_file", type=str, default=None,
                                 help="FASTQ Reads file (-r from bfast).")
        matchparser.add_argument("--bfast_main_indexes", type=int, nargs="*", default=None,
                                 help="Main indexes (-i from bfast).")
        matchparser.add_argument("--bfast_gzipped", action="store_true", default=False,
                                 help="Input file is gzipped (-z from bfast).")
        matchparser.add_argument("--bfast_space", type=int, default=0,
                                 help="0: NT space; 1: color space (-A from bfast) [default: %(default)d]")

        matchparser.set_defaults(func=self.run_bfast_match)

        localalignparser = subparsers.add_parser("localalign", help="localalign help")

        localalignparser.add_argument("--bfast_match_file", type=str, default=None,
                                 help="File from bfast match (-m from bfast).")
        localalignparser.add_argument("--bfast_scoring_matrix_file", type=str, default=None,
                                 help="Scoring matrix file (-x from bfast).")
        localalignparser.set_defaults(func=self.run_bfast_localalign)

        postprocessparser = subparsers.add_parser("postprocess", help="postprocess help")

        postprocessparser.add_argument("--bfast_aligned_file", type=str, default=None,
                                 help="File from bfast localalign (-i from bfast).")
        postprocessparser.add_argument("--bfast_scoring_matrix_file", type=str, default=None,
                                 help="Scoring matrix file (-x from bfast).")
        postprocessparser.add_argument("--bfast_read_group_string", type=str, default=None,
                                       help="Read group string.")

        postprocessparser.set_defaults(func=self.run_bfast_postprocess)

        return parser


    def set_options(self, args):
        """Use args from argparse.parse_args to populate the class.

        """

        self.__init__(args=args)

    # def set_options(self, args):
    #     """Use args from argparse.parse_args to populate the class.

    #     """
        
    #     logger.debug("args = %s" % (args, ))
        
    #     self.__init__(reference_fasta=args.reference_fasta,
    #                   reads_file=args.reads_file,
    #                   main_indexes=args.main_indexes,
    #                   output_file=args.output_file,
    #                   space=args.space,
    #                   threads=args.threads,
    #                   other_params=args.bowtie_other_params)


    def make_match_command(self, args):
        _cmd_string = "%(prog)s match -A %(space)s -t -n %(threads)s -f %(reference_fasta)s -r %(reads_file)s %(other_params)s  > %(output_file)s"

        assert args.bfast_reads_file, "Did not specify input file!"
        assert args.bfast_reference_fasta, "Did not specify FASTA reference file!"
        assert args.bfast_output_file, "Did not specify output file"

        other_params = ""

        if args.bfast_other_params:
            other_params += args.bfast_other_params

        if args.bfast_main_indexes:
            main_indexes = ",".join(map(str, args.bfast_main_indexes))
            other_params += " -i %s" % main_indexes
        
        if args.bfast_gzipped:
            other_params += " -z"

        if args.bfast_temp_dir:
            other_params += " -T %s" % args.bfast_temp_dir

        cmd = _cmd_string % dict(prog=self._cmd, 
                                 reference_fasta=args.bfast_reference_fasta,
                                 reads_file=args.bfast_reads_file,
                                 output_file=args.bfast_output_file,
                                 threads=args.bfast_threads,
                                 space=args.bfast_space,
                                 other_params=other_params)
        
        logger.debug("Command to run: %s" % (cmd, ))

        return cmd

    def make_localalign_command(self, args):
        _cmd_string = "%(prog)s localalign -t -n %(threads)s -f %(reference_fasta)s -m %(match_file)s %(other_params)s  > %(output_file)s"

        assert args.bfast_match_file, "Did not specify input file!"
        assert args.bfast_reference_fasta, "Did not specify FASTA reference file!"
        assert args.bfast_output_file, "Did not specify output file"

        other_params = ""

        if args.bfast_other_params:
            other_params += args.bfast_other_params

        if args.bfast_scoring_matrix_file:
            other_params += " -x %s" % args.bfast_scoring_matrix_file 
        
        cmd = _cmd_string % dict(prog=self._cmd, 
                                 reference_fasta=args.bfast_reference_fasta,
                                 match_file=args.bfast_match_file,
                                 output_file=args.bfast_output_file,
                                 threads=args.bfast_threads,
                                 other_params=other_params)
        
        logger.debug("Command to run: %s" % (cmd, ))

        return cmd

    def make_postprocess_command(self, args):
        _cmd_string = "%(prog)s postprocess -Y 0 -t -n %(threads)s -f %(reference_fasta)s -i %(aligned_file)s %(other_params)s  > %(output_file)s"

        assert args.bfast_aligned_file, "Did not specify input file!"
        assert args.bfast_output_file, "Did not specify output file"

        other_params = ""

        if args.bfast_other_params:
            other_params += args.bfast_other_params

        if args.bfast_scoring_matrix_file:
            other_params += " -x %s" % args.bfast_scoring_matrix_file 
            
        if args.bfast_read_group_string:
            other_params += " -r %s" % args.bfast_read_group_string

        cmd = _cmd_string % dict(prog=self._cmd, 
                                 reference_fasta=args.bfast_reference_fasta,
                                 aligned_file=args.bfast_aligned_file,
                                 output_file=args.bfast_output_file,
                                 threads=args.bfast_threads,
                                 other_params=other_params)
        
        logger.debug("Command to run: %s" % (cmd, ))

        return cmd

    @setup_scratch
    def run_bfast_match(self, args):
        """Run bfast match program from the command line.
        
        Assumes that it is on your PATH.
        
        """

        cmd = self.make_match_command(args)
        logger.debug(cmd)
        utils.safe_run(cmd, shell=False)

    @setup_scratch
    def run_bfast_localalign(self, args):
        """Run bfast localalign program from the command line.
        
        Assumes that it is on your PATH.
        
        """

        cmd = self.make_localalign_command(args)
        utils.safe_run(cmd, shell=False)

    @setup_scratch
    def run_bfast_postprocess(self, args):
        """Run bfast postprocess program from the command line.
        
        Assumes that it is on your PATH.
        
        """
        cmd = self.make_postprocess_command(args)
        utils.safe_run(cmd, shell=False)


def main():    
    bfast = BFASTBase()
    
    usage = "%(prog)s [options] arguments"
    parser = argparse.ArgumentParser(usage=usage)
    parser = bfast.argparse(parser)
    
    args = parser.parse_args()

    args.func(args)

    # bowtierunner.set_options(args)
    
    # bowtierunner.run_bowtie()

if __name__ == "__main__":
    main()
