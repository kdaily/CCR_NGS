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
         r"%s/\1.singles" % config['sickle_params']['output_dir'],
         config['sickle_params'])
def run_sickle(input, output, output_singles, params=None):
    """Run sickle to trim ends of reads based on sequence quality and length.
    
    Also gzips the output to save space.
    
    """
    
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
         r"%s/\1.sai" % config['bwa_aln_params']['output_dir'],
         config['bwa_aln_params'])
def run_bwa_aln(input, output, params=None):
    """Run bwa on individual gzipped fastq files.

    """
    
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
         config['bwa_sampe_params']
         )
def run_bwa_sampe(input, output, output_prefix, sample_name, params=None):
    """Run bwa sampe.


    """
    
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
@transform(run_bwa_sampe, 
           regex(r"(.*).bam"), r"\1.filter.bam",
           config['bamfilter_params'])
def run_filterbam(input, output):
    """
    
    """
    
    params['input'] = input
    params['output'] = output
    
    cmd = "module load %(modules)s\n" % params
    cmd += "samtools view -h -F uUfd -q 1 -b %(input)s > %(output)s" % params
    
    logger.debug("cmd = %s" % (cmd))
    
    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']
    
    job_id = utils.safe_qsub_run(cmd, jobname="filterbam",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@transform(run_filterbam, 
           regex(r"(.*).bam"), r"\1.bam.bai",
           config['bamindex_params'])
def run_indexbam1(input, output, params=None):
    """Run samtools index on bam file.
    
    """
    
    params['input'] = input
    params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "module load %(modules)s\n" % params
    cmd += "samtools index %(input)s %(output)s" % params

    job_id = utils.safe_qsub_run(cmd, jobname="bamindex1",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@follows(run_indexbam1)
@transform(run_filterbam,
           regex(r"(.*).sorted.filter.bam"), r"\1.bam",
           config['cleansam_params'])
def run_cleansam(input, output):
    """Clean up BAM file.
    """
    
    params['input'] = input
    params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']
    
    cmd = "java -Xmx%(maxjheap)s -Djava.io.tmpdir=%(tmp_dir)s -jar %(jar_file)s INPUT=%(input)s OUTPUT=%(output)s" % params

    job_id = utils.safe_qsub_run(cmd, jobname="cleansam",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@transform(run_cleansam,
           regex(r"(.*).bam"), r"\1.bam.bai",
           config['bamindex_params'])
def run_indexbam2(input, output, params=None):
    """Run samtools index on bam file.
    
    """
    
    params['input'] = input
    params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "module load %(modules)s\n" % params
    cmd += "samtools index %(input)s %(output)s" % params

    job_id = utils.safe_qsub_run(cmd, jobname="bamindex2",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@follows(run_indexbam2)
@transform(run_cleansam,
           regex(r"(.*).bam"), r"\1.flagstat.txt",
           config['bamflagstat_params'])
def run_flagstat(input, output, params=None):
    """Run samtools flagstat on bam file.
    
    """
    
    params['input'] = input
    params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "module load %(modules)s\n" % params
    cmd += "samtools flagstat %(input)s > %(output)s" % params

    job_id = utils.safe_qsub_run(cmd, jobname="flagstat",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(20)
@follows(run_cleansam, mkdir(config['picard_markduplicates_params']['output_dir']))
@transform(run_cleansam, 
           regex(r".*/(.*).bam"), r"%s/\1.bam" % config['picard_markduplicates_params']['output_dir'],
           config['picard_markduplicates_params'])
def run_mark_duplicates(input, output, params=None):
    """Set up and run the Picard MarkDuplicates program.
    
    """

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    # # Update input and output from global config object
    # params = config['picard_markduplicates_params']

    params['input'] = input
    params['output'] = output
    params['metrics_file'] = "%s.metrics" % output

    cmd = "module load %(modules)s\n" % params
    cmd += ("java -Xmx%(maxjheap)s -Djava.io.tmpdir=%(tmp_dir)s -jar %(jar_file)s "
            "INPUT=%(input)s OUTPUT=%(output)s METRICS_FILE=%(metrics_file)s " 
            "OPTICAL_DUPLICATE_PIXEL_DISTANCE=%(optical_duplicate_pixel_distance)s" % params)
    
    job_id = utils.safe_qsub_run(cmd, jobname="markdups",
                                 nodes=picard_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))



@jobs_limit(20)
@follows(run_mark_duplicates, mkdir(config['mergebam_params']['output_dir']))
@collate(run_mark_duplicates, 
         regex(r".*/M.(.+).bam"), 
         r"%s/\1.bam" % config['mergebam_params']['output_dir'],
         r"\1",
         config['mergebam_params'])
def run_mergebam(input, output, patient_id=None, params=None):
    """Merge all bam files for a single patient together (this includes matched tumor and normal into one file).
    The separate tumor/normal sample info is encoded in the read groups of the BAM files.
    This is done for future steps (GATK recommends that tumor/normal paired data is run through re-align/re-calibrate step together).
    """
    
        
    # Update input and output from global config object
    params['patient_id'] = patient_id
    params['input'] = " ".join(input)
    params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "module load %(modules)s\n" % params
    cmd += ("java -Xmx%(maxjheap)s -Djava.io.tmpdir=%(tmp_dir)s -jar %(jar_file)s "
            "INPUT=%(input)s OUTPUT=%(output)s SORT_ORDER=%(sort_order)s USE_THREADING=true " % params)

    # cmdline = "--maxjheap=%(maxjheap)s --jar=%(jar_file)s --input %(input)s --output=%(output)s --sort_order=%(sort_order)s MergeSamFiles --use_threading=true" % mergebam_params
    # cmd = "python -m ccrngspy.tasks.Picard %s" % cmdline

    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="mergebam%s" % (params['patient_id']),
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))


@jobs_limit(20)
@follows(run_mk_output_dir,
         mkdir(config['gatk_params']['output_dir']),
         mkdir(config['gatk_realigner_target_creator_params']['output_dir']))
@files(config['gatk_realigner_target_creator_params']['reference_fasta'], 
       config['gatk_realigner_target_creator_params']['output_file'],
       config['gatk_realigner_target_creator_params'])
def run_realign_indel_creator(input, output, params=None):
    """First part of GATK recalibration.

    """

    realign_params['input'] = input
    realign_params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    knowns = ""
    for known in config['gatk_realigner_target_creator_params']['known_files']:
        knowns += "--known %s " % known

    realign_params['knowns'] = knowns
    
    cmd = "module load %(modules)s\n" % params
    cmd += ("java -Xmx%(maxjheap)s -jar %(jar_file)s -nt %(threads)s "
            "-R %(reference_fasta)s -T RealignerTargetCreator "
            "-o %(output)s %(knowns)s" % realign_params)

    job_id = utils.safe_qsub_run(cmd, jobname="realignTC",
                                 nodes=realign_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(20)
@follows(run_indexbam, run_realign_indel_creator, 
         mkdir(config['gatk_indel_realigner_params']['output_dir']))
@transform(run_mark_duplicates,
           regex(r".*/(.*).bam"), 
           r"%s/\1.bam" % config['gatk_indel_realigner_params']['output_dir'],
           config['gatk_indel_realigner_params'])
def run_indel_realigner(input, output, params=None):
    """Second part of GATK recalibration.

    """

    params['input'] = input
    params['output'] = output

    knowns = ""
    for known in config['gatk_realigner_target_creator_params']['known_files']:
        knowns += "-known %s " % known

    params['knowns'] = knowns
    params['target_intervals'] = config['gatk_realigner_target_creator_params']['output_file']
    
    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "module load %(modules)s\n" % params
    cmd += ("java -Xmx%(maxjheap)s -Djava.io.tmpdir=%(tmp_dir)s -jar %(jar_file)s -I %(input)s "
            "-R %(reference_fasta)s -T IndelRealigner -targetIntervals %(target_intervals)s "
            "-o %(output)s %(knowns)s --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4" % params)
    
    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="indelRe",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(20)
@follows(run_indel_realigner, mkdir(config['gatk_base_score_recal_params']['output_dir']))
@transform(run_indel_realigner, regex(r".*/(.*).bam"), r"%s/\1.grp" % config['gatk_base_score_recal_params']['output_dir'])
def run_base_score_recalibrator(input, output, params=None):
    """GATK base score recalibration.

    """
    
    realign_params = config['gatk_base_score_recal_params']
    realign_params['input'] = input
    realign_params['output'] = output

    knowns = ""
    for known in config['gatk_base_score_recal_params']['known_files']:
        knowns += "-knownSites %s " % known

    realign_params['knowns'] = knowns
    
    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "java -Xmx%(maxjheap)s -Djava.io.tmpdir=%(tmp_dir)s -jar %(jar_file)s -T BaseRecalibrator -I %(input)s -R %(reference_fasta)s -o %(output)s %(knowns)s" % realign_params

    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="bsRecal",
                                 nodes=realign_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(20)
@follows(run_base_score_recalibrator)
@transform(run_indel_realigner, 
           regex(r".*/(.*).bam"), 
           r"%s/\1.bam" % config['gatk_base_score_recal_params']['output_dir'], 
           r"%s/\1.grp" % config['gatk_base_score_recal_params']['output_dir'])
def run_write_recalibrated_bam(input, output, bqsr_file, params=None):
    """GATK write BAM file with recalibrated base scores.

    """
    
    recal_params = config['gatk_base_score_recal_params']
    recal_params['input'] = input
    recal_params['output'] = output
    recal_params['bqsr_file'] = bqsr_file
    
    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "java -Xmx%(maxjheap)s -Djava.io.tmpdir=%(tmp_dir)s -jar %(jar_file)s -T PrintReads -I %(input)s -R %(reference_fasta)s -BQSR %(bqsr_file)s -o %(output)s" % recal_params

    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="writerecal",
                                 nodes=recal_params['qsub_nodes'],
                                 # params=recal_params['qsub_params'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))


job_list = [run_mk_output_dir,
            run_bwa_aln,
            run_bwa_sampe,
            run_indexbam1,
            run_cleansam,
            run_flagstat,
            run_mark_duplicates,
            run_mergebam]

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
