#!/usr/bin/env python
"""Run the RNAseq pipeline.

This runs the RNAseq pipeline on a number of files, using qsub to submit to the cluster,
using the version 2 of RUM.

It requires a YAML configuration file with parameters (output directory, etc.)

It also requires a samples file that has at least columns:
    'samplename', 'sample1', 'sample2', 'filename1', and 'filename2'

"""

import argparse
import csv
import sys
import os
import argparse
import subprocess
import time

from ruffus import *
import yaml

from ccrngspy.tasks import Picard
from ccrngspy.tasks import RUM
from ccrngspy.pipeline import rum_helpers
from ccrngspy.pipeline import picard_helpers
from ccrngspy import utils

logger = utils.make_local_logger("Ruffus RNASeq RUM2 QC Logger", level="debug", color=True)

parser = argparse.ArgumentParser(description="Run RUM2 pipeline.")

parser.add_argument('--config_file', dest="config_file", type=str,
                    help="A YAML configuration file for pipeline.")

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

# setup rum specific params
rum_task_params = rum_helpers.make_rum_param_list(samples=samples, config=config, params=None)

#----------------------------------------------
# begin tasks here
#----------------------------------------------

# @follows(mkdir(config['general_params']['log_file_dir']),
#          mkdir(config['fastqc_params']['output_dir']),
#          mkdir(config['rum_params']['output_dir']),
#          mkdir(config['picard_params']['output_dir']))
# def run_setup_dir(input=None, output=None, params=None):
#     """Make high level output directories.

#     """
    
#     pass

# @follows(run_setup_dir)
@follows(mkdir(config['general_params']['log_file_dir']),
         mkdir(config['rum_params']['output_dir']),
         mkdir(config['picard_params']['output_dir']))
def run_mk_output_dir(input=None, output=None, params=None):
    """Make output directories for each sample.

    2012-03-30 Would be better if this was in a mkdir decorator, but not sure how. (KD)
    
    """
    
    if not opts.no_create_output_dir:
        # Make RUM output directory for each sample
        for sample in samples:
            sample_output_dir = os.path.join(config['rum_params']['output_dir'], sample['samplename'])
            try:
                os.mkdir(sample_output_dir)
            except OSError as e:
                logger.debug(e.message)

@follows(run_mk_output_dir)
@files(rum_task_params)
def run_rum(input, output, params=None):
    """Run RUM on paired reads.
    
    """
    
    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
    
    # Add RUM arguments
    rum = RUM.RUMrunner2()
    parser = rum.argparse(parser)
    
    # Update input and output from global config object
    rum_params = config['rum_params']
    rum_params['input'] = input

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['log_file_dir']
    stderr = config['general_params']['log_file_dir']

    ## fastq read files
    rum_params['file1'] = input[0]
    rum_params['file2'] = input[1]
    rum_params['sample'] = params['sample']
    
    cmdline = "--rum_config=%(config)s --rum_run_name=%(sample)s --rum_outdir=%(output_dir)s/%(sample)s --rum_read_files %(file1)s %(file2)s --rum_chunks=%(chunks)s" % rum_params
    args = parser.parse_args(cmdline.split())

    rum.set_options(args)
    
    rum_command = rum.make_command()
    
    # stdout, stderr = utils.safe_run(rum_command, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))

    job_stdout = utils.safe_qsub_run(rum_command, jobname="rum_%s" % params['sample'],
                                     nodes=rum_params['qsub_nodes'], params="-l walltime=168:00:00",
                                     stdout=stdout, stderr=stderr)
    
    logger.debug("stdout = %s" % (job_stdout))

@transform(run_rum, regex(r"(.*).sam"), r"\1.sorted.bam")
def run_sort_sam(input, output, params=None):
    """Set up and run the Picard SortSam program.

    This task works differently than the others; instead of calling the program directly
    by writing out the command line string needed to run it, this runs a python script
    by calling the main function of ccrngspy.tasks.Picard. This is because the Picard code
    is based off of the Galaxy wrapper for Picard, and doesn't work exactly the same as the
    rest.

    2012-03-30 I will consider re-writing it so that it is consistent. (dailykm)
    
    """
    
    # # Let a parser argument handle setting up arguments and options
    # parser = argparse.ArgumentParser()
    
    # # Add Picard arguments
    # picard = Picard.PicardBase()
    # parser = picard.argparse(parser)

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['log_file_dir']
    stderr = config['general_params']['log_file_dir']

    # Update input and output from global config object
    picard_params = config['picard_sortsam_params']

    picard_params['input'] = input
    picard_params['output'] = output

    logger.debug("picard_params = %s" % (picard_params,))
    # Set up using the default arguments, specifying the input and output files since they are required!
    cmdline = "--maxjheap=%(maxjheap)s --jar=%(jar_file)s --input=%(input)s --output=%(output)s --sort_order=%(sort_order)s SortSam" % picard_params

    picard_cmd = "python -m ccrngspy.tasks.Picard %s" % cmdline

    # stdout, stderr = utils.safe_run(picard_cmd, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))
    
    logger.debug("params = %s" % (params, ))
    job_stdout = utils.safe_qsub_run(picard_cmd, jobname="sort_sam",
                                     nodes=picard_params['qsub_nodes'],
                                     stdout=stdout, stderr=stderr)
    
    logger.debug("stdout = %s" % (job_stdout))


@transform(run_sort_sam, regex(r"(.*).bam"), r"\1.bam.bai")
def run_indexbam(input, output, params=None):
    """Run samtools index on bam file.
    
    """

    # params = dict(sample_name=sample_name, index=index, ref_index=ref_index)

    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
        
    # Update input and output from global config object
    bamindex_params = config['bamindex_params']
    bamindex_params['input'] = input
    bamindex_params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['log_file_dir']
    stderr = config['general_params']['log_file_dir']

    cmd = "%(exec)s index %(input)s %(output)s" % bamindex_params

    job_id = utils.safe_qsub_run(cmd, jobname="bamindex",
                                 nodes=bamindex_params['qsub_nodes'], params="-l walltime=168:00:00",
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))



job_list_runfast = [run_mk_output_dir]
job_list_rum = [run_rum]
job_list_rest = [run_sort_sam, run_indexbam]

job_list = job_list_runfast + job_list_rum + job_list_rest

def run_it():
    """Run the pipeline.
    
    Running in three stages to change number of concurrent processes.
    
    1. Run RUM alignment, which can only run on the large memory machine. It would seem that PBS would be smart enough to set the available memory so that concurrent jobs wouldn't happen?
    
    """

    ## Set up directories
    pipeline_run(job_list_runfast, multiprocess=20, logger=logger)

    ## Run RUM
    pipeline_run(job_list_rum, multiprocess=5, logger=logger)

    ## Run sorting
    pipeline_run(job_list_rest, multiprocess=10, logger=logger)

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
        pipeline_printout(sys.stdout, job_list, verbose=3)
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

# run_it()
