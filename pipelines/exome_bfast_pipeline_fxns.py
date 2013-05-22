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
import os
import tempfile

from ruffus import *

from ccrngspy.tasks import BFAST
from ccrngspy import utils

@follows(mkdir(config['general_params']['stdout_log_file_dir']),
         mkdir(config['general_params']['stderr_log_file_dir']))
def run_mk_output_dir(input=None, output=None, params=None):
    """Make stdout and stderr output directories.
    
    """

    pass

@jobs_limit(40)
@follows(run_mk_output_dir, mkdir(config['sickle_params']['output_dir']))
## @collate(sickle_file_list, regex(r".*/(M[NT]..)_(.*)_(R.)_(.*)\.fastq"), r"%s/\1_\2_\4.singles" % config['sickle_params']['output_dir'])
@collate(sickle_file_list, regex(r".*/(M[NT]..)_(.*)_(R.)_(.*)\.fastq.gz"), r"%s/\1_\2_\4.singles" % config['sickle_params']['output_dir'])
def run_sickle(input, output, params=None):
    """Run sickle to trim ends of reads based on sequence quality and length.
    
    """

    # Update input and output from global config object
    params = config['sickle_params']

    input = sorted(input)
    params['input_read1'] = input[0]
    params['input_read2'] = input[1]

    params['input_read1_base'] = os.path.splitext(os.path.basename(input[0]))[0]
    params['input_read2_base'] = os.path.splitext(os.path.basename(input[1]))[0]

    params['output_read1'] = "%(output_dir)s/%(input_read1_base)s" % (params)
    params['output_read2'] = "%(output_dir)s/%(input_read2_base)s" % (params)
    params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    # Modules to load
    
    cmd = ("module load %(modules)s\n"
           "%(exec)s pe -t sanger -l %(length)s -q %(quality)s -f %(input_read1)s -r %(input_read2)s -o %(output_read1)s -p %(output_read2)s -s %(output)s" % params)
    
    logger.debug("cmd = %s" % (cmd,))
    
    job_id = utils.safe_qsub_run(cmd, jobname="sickle",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))
    
@jobs_limit(40)
@follows(run_sickle)
@transform(sickle_gzip_file_list, regex(r".*/(.*).fastq"), r"%s/\1.fastq.gz" % config['sickle_params']['output_dir'])
def run_gzip_sickle(input, output, params=None):
    """Gzip output of sickle files (required for fastqc parsing.)
    
    """
    
    # Update input and output from global config object
    params = config['sickle_params']
    
    params['input'] = input
    params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "gzip %(input)s" % params
    
    logger.debug("cmd = %s" % (cmd,))
    
    job_id = utils.safe_qsub_run(cmd, jobname="gzipsickle",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(40)
@follows(run_gzip_sickle, mkdir(config['merge_paired_reads_params']['output_dir']))
@collate(merge_fastq_file_list, regex(r".*/(.*)_(R[12])_(.*).fastq.gz"), r"%s/\1_\3.fastq.gz" % config['merge_paired_reads_params']['output_dir'])
def run_merge_paired_reads(input, output, params=None):
    """Merge R1 and R2 ends into interleaving file for bfast to work on.
    
    """
    
    # Update input and output from global config object
    params = config['merge_paired_reads_params']
    input = sorted(input)
    params['input_read1'] = input[0]
    params['input_read2'] = input[1]
    params['output'] = output
    
    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir'] 
    
    cmd = "python %(exec)s --gzipped -1 %(input_read1)s -2  %(input_read2)s | gzip -c > %(output)s" % params
    # cmd = "python %(exec)s -1 %(input_read1)s -2  %(input_read2)s | gzip -c > %(output)s" % params
    
    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="mergepairs",
                                 nodes=params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(30)
@follows(run_merge_paired_reads, mkdir(config['bfast_match_params']['output_dir']))
@files(bfast_match_task_params)
def run_bfast_match(input, output, params=None):
    """Run match.
    
    """
    
    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
    
    # Add BFAST arguments
    bfast = BFAST.BFASTBase()
    
    # Use bfast provided by PATH; set using the module functions.
    bfast.set_exec("bfast")
    
    parser = bfast.argparse(parser)
    
    # Update input and output from global config object
    bfast_params = config['bfast_match_params']
    bfast_params['input'] = input
    bfast_params['output'] = output
    
    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']
    
    bfast_params['sample'] = params['sample_name']
    bfast_params['bfast_ref_index'] = params['bfast_ref_index']
    
    cmdline = ("--bfast_temp_dir=%(temp_dir)s --bfast_threads=%(threads)s "
               "--bfast_reference_fasta=%(reference_fasta)s "
               "--bfast_output_file=%(output)s match --bfast_gzipped "
               "--bfast_reads_file=%(input)s --bfast_space=0 "
               "--bfast_main_indexes=%(bfast_ref_index)s" % bfast_params)
    
    args = parser.parse_args(cmdline.split(" "))
    
    bfast_command = ("module load %(modules)s\n"
                     bfast.make_match_command(args) % bfast_params)
    
    logger.debug("cmd = %s" % (bfast_command,))
    
    # stdout, stderr = utils.safe_run(bfast_command, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))
    
    job_id = utils.safe_qsub_run(bfast_command, jobname="bfmatch",
                                 nodes=bfast_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(30)
@follows(run_bfast_match, mkdir(config['bfast_mergematch_params']['output_dir']))
@collate(run_bfast_match, regex(r".*/(.+)_(.+)_(.+)_(.+).bmf"), r"%s/\1_\2_\3.bmf" % config['bfast_mergematch_params']['output_dir'],
           r"\2", r"\3")
def run_merge_bfastmatch(input, output, sample_name=None, index=None):
    """Merge bfast match files using bfast bmfmerge utility with default parameters.

    """
    
    # Update input and output from global config object
    merge_params = config['bfast_mergematch_params']
    merge_params['input'] = " ".join(input)
    merge_params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = ("module load %(modules)s\n"
           "%(exec)s %(input)s > %(output)s" % merge_params)

    job_id = utils.safe_qsub_run(cmd, jobname="mergebfmatch",
                                 nodes=merge_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(30)
@follows(run_merge_bfastmatch, mkdir(config['bfast_localalign_params']['output_dir']))
@transform(run_merge_bfastmatch, regex(r".*/(.*)_(.*)_(.*).bmf"), r"%s/\1_\2_\3.baf" % config['bfast_localalign_params']['output_dir'],
           r"\2", r"\3")
@posttask(touch_file("%s/bfast_localalign_completed.flag" % config['bfast_localalign_params']['output_dir']))
def run_bfast_localalign(input, output, sample_name=None, index=None):
    """Run localalign.
    
    """

    params = dict(sample_name=sample_name, index=index)

    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
    
    # Add BFAST arguments
    bfast = BFAST.BFASTBase()

    # Use bfast provided by PATH; set using the module functions.
    bfast.set_exec("bfast")

    parser = bfast.argparse(parser)
    
    # Update input and output from global config object
    bfast_params = config['bfast_localalign_params']
    bfast_params['input'] = input
    bfast_params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmdline = ("--bfast_threads=%(threads)s "
               "--bfast_reference_fasta=%(reference_fasta)s "
               "--bfast_output_file=%(output)s localalign "
               "--bfast_match_file=%(input)s" % bfast_params)

    args = parser.parse_args(cmdline.split(" "))
    
    bfast_command = ("module load %(modules)s\n"
                     bfast.make_localalign_command(args) % bfast_params)

    # stdout, stderr = utils.safe_run(bfast_command, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))

    job_id = utils.safe_qsub_run(bfast_command, jobname="bfla%s%s" % (params['sample_name'], params['index']),
                                 nodes=bfast_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(30)
@follows(run_bfast_localalign, mkdir(config['bfast_postprocess_params']['output_dir']))
@transform(run_bfast_localalign, regex(r".*/(.*)_(.*)_(.*).baf"), r"%s/\1_\2_\3.sam" % config['bfast_postprocess_params']['output_dir'],
           r"\1", r"\2", r"\3")
@posttask(touch_file("%s/bfast_postprocess_completed.flag" % config['bfast_postprocess_params']['output_dir']))
def run_bfast_postprocess(input, output, flowcell_id=None, sample_name=None, index=None):
    """Run postprocess.
    
    """

    params = dict(sample_name=sample_name, index=index)

    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
    
    # Add BFAST arguments
    bfast = BFAST.BFASTBase()

    # Use bfast provided by PATH; set using the module functions.
    bfast.set_exec("bfast")

    parser = bfast.argparse(parser)
    
    # Update input and output from global config object
    bfast_params = config['bfast_postprocess_params']
    bfast_params['input'] = input
    bfast_params['output'] = output
    read_group_id = read_group_id_dict[sample_name]

    read_group_string = "@RG\tID:%s\tPL:ILLUMINA\tPU:%s\tLB:%s\tSM:%s" % (read_group_id, flowcell_id, config['general_params']['library_name'], sample_name)

    read_group_tempfile = tempfile.NamedTemporaryFile(mode="w", dir=bfast_params['scratch_dir'], delete=True)
    read_group_tempfile.write("%s\n" % read_group_string)
    read_group_tempfile.flush()

    bfast_params['read_group_string'] = read_group_tempfile.name

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']
    
    cmdline = ("--bfast_threads=%(threads)s --bfast_output_file=%(output)s "
               "--bfast_reference_fasta=%(reference_fasta)s postprocess "
               "--bfast_read_group_string=%(read_group_string)s "
               "--bfast_aligned_file=%(input)s" % bfast_params)

    args = parser.parse_args(cmdline.split(" "))
    
    bfast_command = ("module load %(modules)s\n"
                     bfast.make_postprocess_command(args) % bfast_params)

    # stdout, stderr = utils.safe_run(bfast_command, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))

    job_id = utils.safe_qsub_run(bfast_command, jobname="bfpp%s%s" % (params['sample_name'], params['index']),
                                 nodes=bfast_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@jobs_limit(30)
@follows(run_bfast_postprocess, mkdir(config['picard_sortsam_params']['output_dir']))
@transform(run_bfast_postprocess, regex(r".*/(.*)_(.*)_(.*).sam"), r"%s/\1_\2_\3.bam" % config['picard_sortsam_params']['output_dir'],
           r"\2", r"\3")
def run_sort_sam(input, output, sample_name=None, index=None):
    """Run Picard SortSam to convert to sorted bam file.

    """

    params = dict(sample_name=sample_name, index=index)
    
    # Update input and output from global config object
    picard_params = config['picard_sortsam_params']
    picard_params['input'] = input
    picard_params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmdline = "--maxjheap=%(maxjheap)s --jar=%(jar_file)s --input=%(input)s --output=%(output)s --sort_order=%(sort_order)s SortSam " % picard_params

    picard_cmd = "python -m ccrngspy.tasks.Picard %s" % cmdline

    # stdout, stderr = utils.safe_run(picard_cmd, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))
    
    logger.debug("params = %s" % (params, ))
    job_id = utils.safe_qsub_run(picard_cmd, jobname="sortsam",
                                 nodes=picard_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    

    logger.debug("job_id = %s" % (job_id,))
     
@follows(run_sort_sam, mkdir(config['mergebam_params']['output_dir']))
@collate(run_sort_sam, regex(r".*/(.+)_M.(.+)_(.+).bam"), r"%s/\1_\2.bam" % config['mergebam_params']['output_dir'],
         r"\2")
def run_mergebam(input, output, patient_id=None):
    """Merge all bam files for a single patient together (this includes matched tumor and normal into one file).
    The separate tumor/normal sample info is encoded in the read groups of the BAM files.
    This is done for future steps (GATK recommends that tumor/normal paired data is run through re-align/re-calibrate step together).
    """
    
    params = dict(patient_id=patient_id)
        
    # Update input and output from global config object
    mergebam_params = config['mergebam_params']
    mergebam_params['input'] = " ".join(input)
    mergebam_params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmdline = "--maxjheap=%(maxjheap)s --jar=%(jar_file)s --input %(input)s --output=%(output)s --sort_order=%(sort_order)s MergeSamFiles --use_threading=true" % mergebam_params

    cmd = "python -m ccrngspy.tasks.Picard %s" % cmdline
    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="mb%s" % (params['patient_id']),
                                 nodes=mergebam_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

# def run_gsnap(input, output, params=None):
#     """Run gsnap.
    
#     """

#     pass

@transform(run_mergebam, regex(r".*/(.*).bam"), r"%s/\1.bam" % config['picard_markduplicates_params']['output_dir'])
def run_mark_duplicates(input, output, params=None):
    """Set up and run the Picard MarkDuplicates program.

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
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    # Update input and output from global config object
    picard_params = config['picard_markduplicates_params']

    picard_params['input'] = input
    picard_params['output'] = output
    picard_params['metrics_file'] = "%s.metrics" % output

    # Set up using the default arguments, specifying the input and output files since they are required!
    cmdline = "--maxjheap=%(maxjheap)s --jar=%(jar_file)s --input=%(input)s --output=%(output)s MarkDuplicates --metrics_file=%(metrics_file)s --optdupdist=%(optical_duplicate_pixel_distance)s " % picard_params

    picard_cmd = "python -m ccrngspy.tasks.Picard %s" % cmdline

    # stdout, stderr = utils.safe_run(picard_cmd, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))
    
    job_id = utils.safe_qsub_run(picard_cmd, jobname="markdups",
                                 nodes=picard_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))


@transform(run_mark_duplicates, regex(r"(.*).bam"), r"\1.bam.bai")
def run_indexbam(input, output, params=None):
    """Run samtools index on bam file.
    
    """

    # params = dict(sample_name=sample_name, index=index)

    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
        
    # Update input and output from global config object
    bamindex_params = config['bamindex_params']
    bamindex_params['input'] = input
    bamindex_params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "%(exec)s index %(input)s %(output)s" % bamindex_params

    job_id = utils.safe_qsub_run(cmd, jobname="bamindex",
                                 nodes=bamindex_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@follows(run_mk_output_dir, mkdir(config['gatk_realigner_target_creator_params']['output_dir']))
@files(config['gatk_realigner_target_creator_params']['reference_fasta'], config['gatk_realigner_target_creator_params']['output_file'])
def run_realign_indel_creator(input, output, params=None):
    """First part of GATK recalibration.

    """

    realign_params = config['gatk_realigner_target_creator_params']
    realign_params['input'] = input
    realign_params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    knowns = ""
    for known in config['gatk_realigner_target_creator_params']['known_files']:
        knowns += "--known %s " % known

    realign_params['knowns'] = knowns
    
    cmd = "java -Xmx%(maxjheap)s -jar %(jar_file)s -nt %(threads)s -R %(reference_fasta)s -T RealignerTargetCreator -o %(output)s %(knowns)s" % realign_params

    job_id = utils.safe_qsub_run(cmd, jobname="realignTC",
                                 nodes=realign_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@follows(run_indexbam, run_realign_indel_creator, mkdir(config['gatk_indel_realigner_params']['output_dir']))
@transform(run_mark_duplicates, regex(r".*/(.*).bam"), r"%s/\1.bam" % config['gatk_indel_realigner_params']['output_dir'])
def run_indel_realigner(input, output, params=None):
    """Second part of GATK recalibration.

    """

    realign_params = config['gatk_indel_realigner_params']
    realign_params['input'] = input
    realign_params['output'] = output

    knowns = ""
    for known in config['gatk_realigner_target_creator_params']['known_files']:
        knowns += "-known %s " % known

    realign_params['knowns'] = knowns
    realign_params['target_intervals'] = config['gatk_realigner_target_creator_params']['output_file']
    
    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "java -Xmx%(maxjheap)s -Djava.io.tmpdir=%(tmp_dir)s -jar %(jar_file)s -I %(input)s -R %(reference_fasta)s -T IndelRealigner -targetIntervals %(target_intervals)s -o %(output)s %(knowns)s --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4" % realign_params

    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="indelRe",
                                 nodes=realign_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

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

@follows(run_base_score_recalibrator)
@transform(run_indel_realigner, regex(r".*/(.*).bam"), r"%s/\1.bam" % config['gatk_base_score_recal_params']['output_dir'], r"%s/\1.grp" % config['gatk_base_score_recal_params']['output_dir'])
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

@split(run_write_recalibrated_bam, regex(r".*/(.+)_(.+).bam"), r'%s/\1_M[NT]\2.bam' % config['gatk_base_score_recal_params']['output_dir'])
def run_split_bam(input, output, params=None):
    """Split BAM files into separate tumor/normal files based on read group.

    """

    params = config['split_bam_params']
    params['input'] = input
    params['output_dir'] = config['gatk_base_score_recal_params']['output_dir']

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = '%(exec)s %(input)s --output_dir=%(output_dir)s ' % params
    cmd += '--output_file_expr="%(PU)s_%(SM)s.bam"'

    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="splitbam",
                                 nodes=params['qsub_nodes'],
                                 # params=params['qsub_params'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))

@transform(run_split_bam, regex(r"(.*).bam"), r"\1.bam.bai")
def run_index_splitbam(input, output, params=None):
    """Run samtools index on bam file.
    
    """

    # params = dict(sample_name=sample_name, index=index)

    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
        
    # Update input and output from global config object
    bamindex_params = config['bamindex_params']
    bamindex_params['input'] = input
    bamindex_params['output'] = output

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "%(exec)s index %(input)s %(output)s" % bamindex_params

    job_id = utils.safe_qsub_run(cmd, jobname="splitbamindex",
                                 nodes=bamindex_params['qsub_nodes'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))


@jobs_limit(2)
@follows(run_index_splitbam)
@follows(run_split_bam, mkdir(config['mutect_params']['output_dir']))
@collate(run_split_bam, regex(r".*/(.*)_(M[NT])(.*)\.bam"),
         r"%s/\1_\3.txt" % config['mutect_params']['output_dir'],
         # r"%s/\1_\3_coverage.wig" % config['mutect_params']['output_dir'],
         r"%s/\1_\3.vcf" % config['mutect_params']['output_dir'])
def run_mutect(input, output, vcf_file, params=None):
    """Run mutect.

    """

    params = config['mutect_params']

    assert len(input) == 2, "Not the right number of input files (should be 2, received %i)" % len(input)

    # Figure out which file is normal, which is tumor...should be alphabetical?

    input = sorted(input)
    params['input_normal'] = input[0]
    params['input_tumor'] = input[1]

    params['call_stats_file'] = output
    ## params['coverage_file'] = coverage_file
    params['vcf_file'] = vcf_file

    params['log_file'] =  config['mutect_params']['output_dir']

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['stdout_log_file_dir']
    stderr = config['general_params']['stderr_log_file_dir']

    cmd = "%(java_exec)s -Xmx%(maxjheap)s -Djava.io.tmpdir=%(tmp_dir)s -jar %(jar_file)s --analysis_type MuTect --dbsnp %(dbsnp_file)s --cosmic %(cosmic_file)s --input_file:normal %(input_normal)s --input_file:tumor %(input_tumor)s --reference_sequence %(reference_fasta)s --out %(call_stats_file)s --vcf %(vcf_file)s --enable_extended_output" % params

    logger.debug("cmd = %s" % (cmd,))

    job_id = utils.safe_qsub_run(cmd, jobname="mutect",
                                 nodes=params['qsub_nodes'],
                                 params=params['qsub_params'],
                                 stdout=stdout, stderr=stderr)
    
    logger.debug("job_id = %s" % (job_id,))
