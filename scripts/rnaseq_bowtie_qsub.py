#!/usr/bin/env python
"""Run the RNAseq pipeline.

This runs the RNAseq pipeline on a number of files, using qsub to submit to the cluster.

It requires a YAML configuration file with parameters (output directory, etc.)

It also requires a samples file that has at least columns:
    'samplename', 'sample1', 'sample2', 'filename1', and 'filename2'

"""

import argparse
import csv
import sys
import os
import subprocess
import time
import shlex

from ruffus import *
import yaml

from ccrngspy.tasks import FastQC
from ccrngspy.tasks import Picard
from ccrngspy.tasks import Bowtie
from ccrngspy.pipeline import fastqc_helpers
from ccrngspy.pipeline import bowtie_helpers
from ccrngspy.pipeline import picard_helpers
from ccrngspy import utils

logger = utils.make_local_logger("Ruffus RNASeq QC Logger", level="debug", color=True)

parser = argparse.ArgumentParser(description="Run fastqc on files.")

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
with open(os.path.abspath(opts.config_file), 'r') as configfile:
    config = yaml.load(configfile)

# Load the samples tab-separated file
with open(os.path.abspath(opts.sample_file), 'r') as samplefile:
    reader = csv.DictReader(samplefile, delimiter="\t")
    samples = list(reader)

# setup bowtie specific params
bowtie_task_params = bowtie_helpers.make_bowtie_param_list(samples=samples, config=config, params=None)


#----------------------------------------------
# begin tasks here
#----------------------------------------------
@follows(mkdir(config['general_params']['log_file_dir']),
         mkdir(config['bowtie_params']['output_dir']),
         mkdir(config['picard_params']['output_dir']))
def run_mk_output_dir(input=None, output=None, params=None):
    """Make output directories for each sample.

    2012-03-30 Would be better if this was in a mkdir decorator, but not sure how. (KD)
    
    """
    
    if not opts.no_create_output_dir:
        # Make Bowtie output directory for each sample
        for sample in samples:
            sample_output_dir = os.path.join(config['bowtie_params']['output_dir'], sample['samplename'])

            try:
                os.mkdir(sample_output_dir)
            except OSError:
                logger.debug("Directory already exists, continuing.")

@follows(run_mk_output_dir)
@files(bowtie_task_params)
def run_bowtie(input, output, params=None):
    """Run Bowtie on paired reads.
    
    """
    
    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
    
    # Add Bowtie arguments
    bowtie = Bowtie.BowtieRunner()
    parser = bowtie.argparse(parser)
    
    # Update input and output from global config object
    bowtie_params = config['bowtie_params']
    bowtie_params['input'] = input

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['log_file_dir']
    stderr = config['general_params']['log_file_dir']

    ## fastq read files
    bowtie_params['file1'] = input[0]
    bowtie_params['file2'] = input[1]
    bowtie_params['sample'] = params['sample']
    bowtie_params['output'] = output

    logger.debug('bowtie_params = %s' % (bowtie_params, ))
    
    cmdline = '--bowtie_index=%(bowtie_index)s -1 %(file1)s -2 %(file2)s -o %(output)s --threads=%(threads)s --other_params="%(other_params)s"' % bowtie_params
    # bowtie_cmd = "python -m ccrngspy.tasks.Bowtie %s" % cmdline

    args = parser.parse_args(shlex.split(cmdline))
    logger.debug("cmdline = %s" % (shlex.split(cmdline), ))


    bowtie.set_options(args)
    
    bowtie_command = bowtie.make_command()
    
    logger.debug("cmd = %s" % (bowtie_command, ))
    logger.debug("params = %s" % (params, ))
    
    # job_stdout, job_stderr = utils.safe_qsub_run(bowtie_command, jobname="bowtie_%s" % params['sample'],
    #                                              nodes=bowtie_params['qsub_nodes'],
    #                                              params="-v np=%(threads)s" % bowtie_params,
    #                                              stdout=stdout, stderr=stderr)
    # logger.debug("stdout = %s, stderr = %s" % (job_stdout, job_stderr))

    job_stdout = utils.safe_qsub_run(bowtie_command, jobname="bowtie_%s" % params['sample'],
                                     nodes=bowtie_params['qsub_nodes'],
                                     params="-v np=%(threads)s" % bowtie_params,
                                     stdout=stdout, stderr=stderr)
    
    logger.debug("stdout = %s" % (job_stdout))

@transform(run_bowtie, regex(r"(.*).sam"), r"\1.sorted.sam")
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

@transform(run_sort_sam, regex(r".*/(.*)/bowtie2.sorted.sam"), r"%s/\1.tsv" % config['picard_params']['output_dir'], r"\1")
def run_collect_rnaseq_metrics(input, output, sample):
    """Set up and run the Picard CollectRnaSeqMetrics program.

    This task works differently than the others; instead of calling the program directly
    by writing out the command line string needed to run it, this runs a python script
    by calling the main function of ccrngspy.tasks.Picard. This is because the Picard code
    is based off of the Galaxy wrapper for Picard, and doesn't work exactly the same as the
    rest.

    2012-03-30 I will consider re-writing it so that it is consistent. (dailykm)

    """
    
    # # Let a parser argument handle setting up arguments and options
    # parser = argparse.ArgumentParser()
    
    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['log_file_dir']
    stderr = config['general_params']['log_file_dir']

    # Update input and output from global config object
    picard_params = config['picard_params']
    picard_params['input'] = input
    picard_params['output'] = output
    
    # Set up using the default arguments, specifying the input and output files since they are required!
    cmdline = "--maxjheap=%(maxjheap)s --jar=%(jar_file)s --input=%(input)s --output=%(output)s --ref_flat=%(ref_flat)s --ref_file=%(ref_file)s CollectRnaSeqMetrics --minimum_length=%(minimum_length)s --chart_output=%(chart_output)s --metric_accumulation_level=%(metric_accumulation_level)s --stop_after=%(stop_after)s" % picard_params

    # args = parser.parse_args(cmdline.split())
    
    # # Run the function for collecting RNASeq metrics
    # args.func(args)
    
    picard_cmd = "python -m ccrngspy.tasks.Picard %s" % cmdline

    # stdout, stderr = utils.safe_run(picard_cmd, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))

    job_stdout = utils.safe_qsub_run(picard_cmd, jobname="rnaseqmet_%s" % sample,
                                     nodes=picard_params['qsub_nodes'],
                                     stdout=stdout, stderr=stderr)
    
    logger.debug("stdout = %s" % (job_stdout))

@merge(run_collect_rnaseq_metrics, os.path.join(config["picard_params"]["output_dir"], "CollectRNASeqMetrics.tsv"))
def run_merge_rnaseq_metrics(input_files, summary_file):
    """Merge the outputs of collectrnaseqmetrics into one tab-separated file.

    """

    metrics = []
  
    for fn in input_files:
        ## Only take the first non-comment line; there is some histogram data after it!
   
        metrics.extend(picard_helpers.parse_picard_rnaseq_metrics(fn)[0:1])

    fieldnames = metrics[0].keys()
    
    with open(summary_file, 'w') as fou:
        dw = csv.DictWriter(fou, delimiter='\t', fieldnames=fieldnames)
        dw.writeheader()
        dw.writerows(metrics)

job_list_runfast = [run_mk_output_dir]
job_list_rest = [run_bowtie, run_sort_sam, run_collect_rnaseq_metrics, run_merge_rnaseq_metrics]

def run_it():
    """Run the pipeline.
    
    Running in three stages to change number of concurrent processes.
    
    1. Run Bowtie2 alignment.
    2. Run CollectRNASeqMetrics etc, which are memory intensive.
    
    """

    ## Set up directories
    pipeline_run(job_list_runfast, multiprocess=20, logger=logger)
    
    ## Run Bowtie, sorting, collecting RNASeq metrics
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
