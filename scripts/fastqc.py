#!/usr/bin/env python
"""Run the fastqc pipeline.

This runs FastQC on a number of files.

It requires a YAML configuration file with parameters for FastQC (output directory, etc.)
It also requires a samples file that has at least a column named 'sample' and 'filename'.

"""

import argparse
import csv
import sys
import os
import argparse
import subprocess
import time
import tempfile
from collections import defaultdict

from ruffus import *
import yaml

from ccrngspy.tasks import FastQC
from ccrngspy import utils

from ccrngspy.pipeline import fastqc_helpers
from ccrngspy.pipeline import dummy_helpers

logger = utils.make_local_logger("FastQC logging", level="debug", color="green")

parser = argparse.ArgumentParser(description="Run fastqc on files.")

parser.add_argument("--print_only", dest="print_only", action="store_true", default=False,
                    help="Don't run the pipeline, just print what will be run.")

parser.add_argument('--verbose', type=int, default=3,
                    help="Verbosity when using print only mode.")

parser.add_argument('--config_file', dest="config_file", type=str,
                    help="A YAML configuration file for pipeline.")

parser.add_argument('--sample_file', dest="sample_file", type=str,
                    help="A YAML configuration file for pipeline.")

# add options for the fastqc task
parser = FastQC.FastQC().argparse(parser)

# Parse the options
opts = parser.parse_args()

# Load the bootstrap config file
with open(opts.config_file, 'r') as configfile:
    config = yaml.load(configfile)

# Load the samples tab-separated file
with open(opts.sample_file, 'r') as samplefile:
    reader = csv.DictReader(samplefile, delimiter="\t")
    samples = list(reader)


fastqc_input_files = fastqc_helpers.make_fastqc_input_list(samples=samples, config=config)

# Subset for testing
fastqc_input_files = fastqc_input_files[:36]

#----------------------------------------------
# begin tasks here
#----------------------------------------------

@collate(fastqc_input_files,
         regex(r".*/(M[NT]..)_(.*)_(R.)_.*\.fastq.gz"),
         r"%s/\1_\2_\3_fastqc.zip" % config['fastqc_params']['output_dir'])
def run_fastqc(input, output, params=None):
    """Set up and run the fastqc program.
    
    """

    params = config['fastqc_params']

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']
    
    try:
        tmp = params['casava']
    except KeyError as e:
        logger.info("Casava parameter not specified, assuming false.")
        params['casava'] = False
        
        
    fastqc_task = FastQC.FastQC(input_files=input, output_directory=params['output_dir'],
                                casava=params['casava'], threads=params['threads'])

    cmd = fastqc_task.make_command()

    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="fastqc",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))


job_list =  [run_fastqc]

def run_it():
    """Run the pipeline.
        
    """

    ## Run fastqc
    pipeline_run(job_list, multiprocess=20, logger=logger)


def _keep_alive():
    """Do something easy so that any task-killing program thinks that I am still alive!

    """

    for x in xrange(10000):
        foo = max(y for y in xrange(1000))
    
    time.sleep(10)


## In order to fool any task-killing software, we'll run the pipeline in one thread (that may be dormant for a while)
## and some simple keep-alive task in another.
import multiprocessing

if opts.print_only:
    pipeline_printout(sys.stdout, job_list, verbose=opts.verbose)
else:
    print "Starting the main program"
    thread = multiprocessing.Process(target=run_it)

    print "Launching Pipeline"
    thread.start()
    print "Pipeline has been launched"
    
    # is_alive() is False when the thread ends.
    while thread.is_alive():
        
        # Here you would have the code to make the app
        # look "alive", a progress bar, or maybe just
        # keep on working as usual.
        
        _keep_alive()

        # print "The pipeline is still running"
        # sys.stdout.flush()

        # Wait a little bit, or until the thread ends,
        # whatever's shorter.
        thread.join(0.3)

    print "Pipeline done!"


