#!/usr/bin/env python
"""Run a BWA exome pipeline, using qsub to submit to the cluster.

The FASTQ files must be preprocessed to be in the format required for paired end sequencing.

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
from ccrngspy.pipeline import bwa_helpers
from ccrngspy import utils

logger = utils.make_local_logger("Ruffus BWA QC Logger", level="debug", color=True)

parser = argparse.ArgumentParser(description="Run BWA pipeline on files.")

parser.add_argument('--config_file', dest="config_file", type=str,
                    help="A YAML configuration file for pipeline.")

parser.add_argument('--verbose', type=int, default=3,
                    help="Verbosity when using print only mode.")

parser.add_argument('--sample_file', dest="sample_file", type=str,
                    help="A tab separated file about the samples to run.")

parser.add_argument("--print_only", dest="print_only", action="store_true", default=False,
                    help="Don't run the pipeline, just print what will be run.")

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
if config['general_params']['samples_to_run']:
    runthesesamples = config['general_params']['samples_to_run']
    samples = [x for x in samples if x['sample_name'] in runthesesamples]
    
# setup inital run params
sickle_file_list = bwa_helpers.make_sickle_file_list(samples=samples, config=config, params=None)
sickle_outputs = bwa_helpers.make_bwa_files(samples, config)

print sickle_file_list

@follows(mkdir(config['general_params']['stdout_log_file_dir']),
         mkdir(config['general_params']['stderr_log_file_dir']))
def run_mk_output_dir(input=None, output=None, params=None):
    """Make stdout and stderr output directories.
    
    """

    pass

# /data/dailykm/CSAS13444/fastq/D1JF9ACXX/merged/MN19_CAGATC_L002_R2.fastq.gz

@jobs_limit(40)
@follows(run_mk_output_dir, mkdir(config['sickle_params']['output_dir']))
@collate(sickle_file_list, regex(r".*/(M[NT]..)_(R.)\.fastq.gz"),
         [r"%s/\1_R1.fastq.gz" % config['sickle_params']['output_dir'],
          r"%s/\1_R2.fastq.gz" % config['sickle_params']['output_dir']],
         r"%s/\1.singles" % config['sickle_params']['output_dir'])
def run_sickle(input, output, output_singles, params=None):
    """Run sickle to trim ends of reads based on sequence quality and length.
    
    Also gzips the output to save space.
    
    """

    # Update input and output from global config object
    params = config['sickle_params']

    input = sorted(input)
    params['input_read1'] = input[0]
    params['input_read2'] = input[1]

    params['input_read1_base'] = os.path.splitext(os.path.basename(input[0]))[0]
    params['input_read2_base'] = os.path.splitext(os.path.basename(input[1]))[0]

    params['output_read1'] = output[0].split(".gz")[0] # "%(output_dir)s/%(input_read1_base)s" % (params)
    params['output_read2'] = output[1].split(".gz")[0] # "%(output_dir)s/%(input_read1_base)s" % (params)
    params['output_singles'] = output_singles

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = ("module load %(modules)s\n"
           "%(exec)s pe -t sanger -l %(length)s -q %(quality)s -f %(input_read1)s -r %(input_read2)s -o %(output_read1)s -p %(output_read2)s -s %(output_singles)s\n"
           "gzip %(output_read1)s &\n" 
           "gzip %(output_read2)s &\n" 
           "wait" % params)
    
    logger.debug("cmd = %s" % (cmd,))
    
    job_id = utils.safe_qsub_run(cmd, jobname="sickle",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@follows(run_sickle,
         mkdir(config['bwa_aln_params']['output_dir']))
@collate(sickle_outputs,
         regex(r".*/(.*).fastq.gz"), 
         r"%s/\1.sai" % config['bwa_aln_params']['output_dir'])
def run_bwa_aln(input, output, params=None):
    """Run bwa on individual gzipped fastq files.

    """
    
    # Update input and output from global config object
    params = config['bwa_aln_params']
    params['input'] = " ".join(input)
    params['output'] = output
    
    cmd = "module load %(modules)s\n" % params
    cmd += "bwa aln -t %(threads)s %(reference_fasta)s %(input)s > %(output)s" % params

    logger.debug("cmd = %s" % (cmd))

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']
    
    job_id = utils.safe_qsub_run(cmd, jobname="bwa_aln",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@follows(run_bwa_aln,
         mkdir(config['bwa_sampe_params']['output_dir']))
@collate(run_bwa_aln,
         regex(r".*/(.*)_(.*).sai"), 
         add_inputs([r"%s/\1_R1.fastq.gz" % config['sickle_params']['output_dir'],
                     r"%s/\1_R2.fastq.gz" % config['sickle_params']['output_dir']]),
         r"%s/\1.sorted.bam" % config['bwa_sampe_params']['output_dir'],
         r"\1.sorted",
         r"\1",
         )
def run_bwa_sampe(input, output, output_prefix, sample_name, params=None):
    """Run bwa sampe.


    """
    
    # Update input and output from global config object
    params = config['bwa_sampe_params']
    
    params['sai_R1'] = input[0][0]
    params['sai_R2'] = input[1][0]

    # For some reason, -f param didn't work with samtools sort
    # So, need to use prefix verison (without bam suffix)
    params['output'] = output_prefix
    
    (params['fastq_R1'], params['fastq_R2']) = input[0][1]
        
    params['read_group_string'] = bwa_helpers.sample_name_to_read_group_string(sample_name=sample_name, samples=samples, config=config)
    cmd = "module load %(modules)s\n" % params
    cmd = ('bwa sampe -r "%(read_group_string)s" %(reference_fasta)s %(sai_R1)s %(sai_R2)s %(fastq_R1)s %(fastq_R2)s | '
           'samtools view -Su - | '
           'samtools sort - -@ %(threads)s -f %(output)s' % params)
    
    logger.debug("cmd = %s" % (cmd))
    
    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']
    
    job_id = utils.safe_qsub_run(cmd, jobname="bwa_sampe",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@follows(run_bwa_sampe,
         mkdir(config['bamfilter_params']['output_dir']))
@transform(run_bwa_sampe, regex(r"(.*).bam"), r"\1.filter.bam")
def run_filterbam(input, output):
    """
    
    """

    # Update input and output from global config object
    params = config['bamfilter_params']

    params['input'] = input
    params['output'] = output
    
    cmd = "module load %(modules)s\n" % params
    cmd = "samtools view -h -F uUfd -q 1 -b %(input)s > %(output)s" % params
    
    logger.debug("cmd = %s" % (cmd))
    
    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']
    
    job_id = utils.safe_qsub_run(cmd, jobname="filterbam",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@transform(run_filterbam, regex(r"(.*).bam"), r"\1.bam.bai")
def run_indexbam1(input, output, params=None):
    """Run samtools index on bam file.
    
    """


    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
        
    # Update input and output from global config object
    params = config['params']
    params['input'] = input
    params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "samtools index %(input)s %(output)s" % params

    job_id = utils.safe_qsub_run(cmd, jobname="bamindex",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

job_list = [run_mk_output_dir, run_bwa_aln, run_bwa_sampe, run_indexbam1]

def run_it():
    """Run the pipeline.
        
    """
    
    pipeline_run(job_list, multiprocess=40, logger=logger, gnu_make_maximal_rebuild_mode=False)

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
