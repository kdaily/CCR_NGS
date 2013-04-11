"""Helper files for pipeline bfast steps.

"""

import os
from collections import defaultdict

from ccrngspy import utils

logger = utils.make_local_logger("BFAST helper logging", level="debug", color="green")

def make_bfast_readgroup_lookup(samples, config, params=None):
    """Helper function to make a lookup dictionary to find which read group ID to set.

    This is currently used during the bfast postprocess step, when the read group is added.

    The current setup now is that the ID is based on the sample, so when tumor/normal pairs are merged,
    the read groups are separated by sample.

    The read_group_id must be a column in the samples file, called read_group_id.
    
    """

    read_group_id_dict = {}
    
    for sample in samples:
        read_group_id_dict[sample['sample_name']] = sample['read_group_id']
        
    return read_group_id_dict
    
    
def make_sickle_param_list(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    Needs to be a list of [input, output, params]; for the fastqc script.

    The output is determined by the flowcell id, sample name, and reference fasta index.

    while the params are taken from the global opts variable
    (and possibly from the YAML config file).
    
    """

    final_list = []

    fastq_dir = config['general_params']['fastq_input_dir']
    output_dir = config['sickle_params']['output_dir']
    
    for sample in samples:
        params = dict(sample_name=sample['sample_name'], index=sample['index'])
            
        # output_file = "%(flowcell_id)s_%(sample_name)s_%(tag)s_%(lane)s_%(index)s.fastq.gz" % sample
            
        tmp1 = [os.path.join(fastq_dir, sample['filename1']),
                os.path.join(output_dir, sample['filename1']),
                params]

        tmp2 = [os.path.join(fastq_dir, sample['filename2']),
                os.path.join(output_dir, sample['filename2']),
                params]

        final_list.extend([tmp1, tmp2])

    return final_list

def make_sickle_file_list(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    Needs to be a list of [input]; for the sickle script.
    
    """

    final_list = []

    fastq_dir = config['general_params']['fastq_input_dir']
    
    for sample in samples:
        final_list.extend([os.path.join(fastq_dir, sample['filename1']), os.path.join(fastq_dir, sample['filename2'])])

    return final_list

def make_merge_fastq_file_list(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    Needs to be a list of [input]; for the sickle script.
    
    """

    final_list = []

    fastq_dir = config['sickle_params']['output_dir']
    
    for sample in samples:
        final_list.extend([os.path.join(fastq_dir, sample['filename1']), os.path.join(fastq_dir, sample['filename2'])])

    return final_list

def make_gzip_fastq_file_list(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    Needs to be a list of [input]; for gzipping.
    
    """

    final_list = []

    fastq_dir = config['sickle_params']['output_dir']
    
    for sample in samples:
        final_list.extend([os.path.join(fastq_dir, os.path.splitext(sample['filename1'])[0]), os.path.join(fastq_dir, os.path.splitext(sample['filename2'])[0])])

    return final_list

def make_bfast_match_param_list(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    Needs to be a list of [input, output, params]; for the fastqc script.

    The output is determined by the flowcell id, sample name, and reference fasta index.

    while the params are taken from the global opts variable
    (and possibly from the YAML config file).
    
    """

    final_list = []

    fastq_dir = config['merge_paired_reads_params']['output_dir']
    output_dir = config['bfast_match_params']['output_dir']
    
    for sample in samples:
        for bfast_ref_index in config['bfast_match_params']['bfast_reference_fasta_indexes']:
            
            params = dict(sample_name=sample['sample_name'], index=sample['index'], bfast_ref_index=bfast_ref_index)

            sample['bfast_ref_index'] = bfast_ref_index
            
            input_file = "%(sample_name)s_%(tag)s_%(lane)s_%(index)s.fastq.gz" % sample
            output_file = "%(flowcell_id)s_%(sample_name)s_%(index)s_%(bfast_ref_index)s.bmf" % sample
            
            tmp = [os.path.join(fastq_dir, input_file),
                   os.path.join(output_dir, output_file),
                   params]
        
            final_list.append(tmp)

    return final_list
