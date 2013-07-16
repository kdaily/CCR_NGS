#!/usr/bin/env python
"""Run the Breakdancer pipeline.

This runs the Breakdancer pipeline on a number of files, using qsub to submit to the cluster.

Requires a directory of BAM files.

It requires a YAML configuration file with parameters (output directory, etc.)

It also requires a samples file that has at least columns:
'source, 'samplename', 'filename'

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

# from ccrngspy.tasks import BFAST
# from ccrngspy.pipeline import bfast_helpers
from ccrngspy import utils

logger = utils.make_local_logger("Ruffus Breakdancer QC Logger", level="debug", color=True)

parser = argparse.ArgumentParser(description="Run breakdancer on BAM files.")

parser.add_argument('--config_file', dest="config_file", type=str,
                    help="A YAML configuration file for pipeline.")

parser.add_argument('--verbose', type=int, default=3,
                    help="Verbosity when using print only mode.")

parser.add_argument('--sample_file', dest="sample_file", type=str,
                    help="A tab separated file about the samples to run.")

parser.add_argument("--print_only", dest="print_only", action="store_true", default=False,
                    help="Don't run the pipeline, just print what will be run.")

parser.add_argument("--no_output_dir", dest="no_create_output_dir", action="store_true", default=False,
                    help="Don't recreate the output dirs.")

# add options for the fastqc task
# parser = FastQC.FastQC().argparse(parser)

# Parse the options
opts = parser.parse_args()

# Load the bootstrap config file
with open(opts.config_file, 'r') as configfile:
    config = yaml.load(configfile)

# Load the samples tab-separated file
with open(opts.sample_file, 'r') as samplefile:
    reader = csv.DictReader(samplefile, delimiter="\t")
    samples = list(reader)


## Only this sample should run
runthesesamples = config['general_params']['samples_to_run'] 
samples = [x for x in samples if x['sample_name'] in runthesesamples]

print samples
    
# setup inital run params
bam_file_list = [os.path.join(config['general_params']['bam_input_dir'], x['filename']) for x in samples]

@follows(mkdir(config['general_params']['stdout_log_file_dir']),
         mkdir(config['general_params']['stderr_log_file_dir']))
def run_mk_output_dir(input=None, output=None, params=None):
    """Make stdout and stderr output directories.
    
    """

    pass

@jobs_limit(40)
@follows(run_mk_output_dir, 
         mkdir(config['bam2cfg_params']['output_dir']))
@transform(bam_file_list, 
           regex(r".*/(.*).bam"), 
           r"%s/\1.cfg" % config['bam2cfg_params']['output_dir'],
           config['bam2cfg_params']
           )
def run_bam2cfg(input, output, params=None):
    """Create breakdancer config files.
    
    """

    # Update input and output from global config object
    params['input'] = input
    params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "module load %(modules)s\n" % params
    cmd += "%(exec)s -g -h %(input)s > %(output)s" % params
    
    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="bam2cfg",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(40)
@follows(run_bam2cfg, mkdir(config['breakdancer_params']['output_dir']))
@transform(run_bam2cfg, 
           regex(r".*/(.*).cfg"), 
           r"%s/\1.txt" % config['breakdancer_params']['output_dir'], 
           r"%s/\1.bed" % config['breakdancer_params']['output_dir'],
           r"%s/\1.fastq" % config['breakdancer_params']['output_dir'],
           config['breakdancer_params'])
def run_breakdancer(input, output, bedfile, fastqfile, params=None):
    """Run breakdancer, only to find intra-chromosomal structural variations (-t param).
    
    """

    params['input'] = input
    params['output'] = output
    params['bedfile'] = bedfile
    params['fastqfile'] = fastqfile

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "module load %(modules)s\n" % params
    cmd += "%(exec)s -g %(bedfile)s -d %(fastqfile)s %(transchrom)s -r %(min_reads)s %(input)s > %(output)s" % params
    
    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="breakdancer",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

job_list = [run_bam2cfg, run_breakdancer]

def run_it():
    """Run the pipeline.
        
    """
    
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
