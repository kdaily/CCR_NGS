#!/usr/bin/env python
"""Run the BFAST exome pipeline.

This runs the BFAST pipeline on a number of files, using qsub to submit to the cluster.

The FASTQ files must be preprocessed to be in the format required for paired end sequencing.
The paired reads are in the same file, consecutively (i.e. seq1_R1, seq1_R2, seq2_R1, seq2_R2, etc.)

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

from ccrngspy.tasks import BFAST
from ccrngspy.pipeline import bfast_helpers
from ccrngspy import utils

logger = utils.make_local_logger("Ruffus BFAST QC Logger", level="debug", color=True)

parser = argparse.ArgumentParser(description="Run bfast on files.")

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
runthesesamples = ['MN24', 'MT24', 'MN23', 'MT23', 'MN22', 'MT22', 'MN21', 'MT21', 'MN20', 'MT20', 'MN19', 'MT19']
samples = [x for x in samples if x['sample_name'] in runthesesamples]
    
# setup inital run params
sickle_file_list = bfast_helpers.make_sickle_file_list(samples=samples, config=config, params=None)

# files to gzip
sickle_gzip_file_list = bfast_helpers.make_gzip_fastq_file_list(samples=samples, config=config, params=None)

# setup merge fastq for bfast params
merge_fastq_file_list = bfast_helpers.make_merge_fastq_file_list(samples=samples, config=config, params=None)

# Uses sickle's output dir to setup input file location
bfast_match_task_params = bfast_helpers.make_bfast_match_param_list(samples=samples, config=config, params=None)

# set up read group id dict for lookup
read_group_id_dict = bfast_helpers.make_bfast_readgroup_lookup(samples=samples, config=config, params=None)

from ccrngspy.ruffus_pipelines.exome_bfast_pipeline_fxns import *
    
job_list_runfast = [run_mk_output_dir, run_sickle, run_gzip_sickle, run_merge_paired_reads]
job_list_bfast = [run_bfast_match, run_merge_bfastmatch, run_bfast_localalign, run_bfast_postprocess]
job_list_postprocess = [run_sort_sam, run_mergebam, run_mark_duplicates, run_indexbam]
job_list_GATK = [run_realign_indel_creator, run_indel_realigner, run_base_score_recalibrator, run_write_recalibrated_bam]
job_list_mutations = [run_split_bam, run_mutect]

job_list = job_list_runfast + job_list_bfast + job_list_postprocess + job_list_GATK + job_list_mutations

def run_it():
    """Run the pipeline.
        
    """
    
    pipeline_run(job_list, multiprocess=20, logger=logger)
    
    # ## Set up directories
    # pipeline_run(job_list_runfast, multiprocess=20, logger=logger)

    # ## Run BFAST
    # pipeline_run(job_list_bfast, multiprocess=20, logger=logger)

    # ## Run sorting, collecting RNASeq metrics
    # pipeline_run(job_list_rest, multiprocess=2, logger=logger)

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
