"""Wrapper for Tophat.

Kenny Daily, 2012

"""

import subprocess
import argparse
import os
import Task

from ccrngspy import utils

logger = utils.make_local_logger("Tophat logging", level="debug", color="green")

class TophatRunner(Task.Task):
    """Container for Tophat tasks.
    
    """
    
    _cmd = "/usr/local/tophat-2.0.4/bin/tophat2"
    
    def __init__(self, index=None, mate_file_one=None, mate_file_two=None, output_file=None, threads=1, other_params=None):

        self.index = index
        self.mate_file_one = mate_file_one
        self.mate_file_two = mate_file_two
        self.output_file = output_file

        self.threads = threads
        self.other_params = other_params
        
    def argparse(self, parser):
        """Add Tophat option group to an OptionParser.
        
        """

        group = parser.add_argument_group("Tophat Options")
        group.add_argument("--bowtie_index", type=str, default=None,
                           help="Bowtie index.")
        group.add_argument("-1", "--mate_file_one", dest="mate_file_one", type=str, default=None,
                           help="First mated read file.")
        group.add_argument("-2", "--mate_file_two", dest="mate_file_two", type=str, default=None,
                           help="Second mated read file.")
        group.add_argument("-o", "--tophat_output_file", type=str, default=None,
                           help="Tophat output file.")
        group.add_argument("--threads", type=int, default=1,
                           help="Specifies the number of threads for Tophat.")
        group.add_argument("--other_params", dest="tophat_other_params", type=str, default=None,
                           help="Other tophat params.")
        return parser
    
    def set_options(self, args):
        """Use args from argparse.parse_args to populate the class.

        """
        
        logger.debug("args = %s" % (args, ))
        
        self.__init__(index=args.bowtie_index,
                      mate_file_one=args.mate_file_one,
                      mate_file_two=args.mate_file_two,
                      output_file=args.tophat_output_file,
                      threads=args.threads,
                      other_params=args.tophat_other_params)


    def make_command(self):
        _cmd_string = "%(prog)s -x %(index)s %(mate_file_one)s %(mate_file_two)s -S %(output_file)s -p %(threads)s %(other_params)s"

        assert (self.mate_file_one and self.mate_file_two), "Did not specify input files!"
        assert self.index, "Did not specify index file!"
        assert self.output_file, "Did not specify output file"

        cmd = _cmd_string % dict(prog=self._cmd, 
                                 index=self.index,
                                 mate_file_one=self.mate_file_one,
                                 mate_file_two=self.mate_file_two,
                                 output_file=self.output_file,
                                 threads=self.threads,
                                 other_params=self.other_params)
        
        logger.debug("Command to run: %s" % (cmd, ))

        return cmd
    
    def run_tophat(self):
        """Run the tophat program from the command line.
        
        Assumes that it is on your PATH.
        
        """

        cmd = self.make_command()
        utils.safe_run(cmd, shell=False)

def main():    
    _test()

def _test():

    tophatrunner = TophatRunner()

    usage = "%(prog)s [options] input_files"
    parser = argparse.ArgumentParser(usage=usage)
    parser = tophatrunner.argparse(parser)
    
    args = parser.parse_args()

    tophatrunner.set_options(args)
    
    tophatrunner.run_tophat()

if __name__ == "__main__":
    main()
