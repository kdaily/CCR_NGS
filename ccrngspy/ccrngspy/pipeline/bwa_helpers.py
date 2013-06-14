"""Helper files for pipeline bwa steps.

"""

import os
from collections import defaultdict

from ccrngspy import utils

logger = utils.make_local_logger("BWA helper logging", level="debug", color="green")

def make_sickle_file_list(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    Needs to be a list of [input]; for the sickle script.
    
    """

    final_list = []

    fastq_dir = config['general_params']['fastq_input_dir']
    
    for sample in samples:
        final_list.extend([os.path.join(fastq_dir, sample['filename1']),
                           os.path.join(fastq_dir, sample['filename2'])])

    return final_list

def make_bwa_files(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    """

    final_list = []

    fastq_dir = config['sickle_params']['output_dir']
    
    for sample in samples:
        
        input_file1 = "%(sample_name)s_R1.fastq.gz" % sample
        input_file2 = "%(sample_name)s_R2.fastq.gz" % sample
        final_list.append(os.path.join(fastq_dir, input_file1))
        final_list.append(os.path.join(fastq_dir, input_file2))

    return final_list

def sample_name_to_read_group_string(sample_name, samples, config):
    read_group_string_template = "@RG\tID:%(sample_name)s\tPL:ILLUMINA\tPU:%(flowcell_id)s\tLB:%(sample_name)s\tSM:%(sample_name)s\tDS:%(library_name)s\tCN:%(sequencing_center)s"

    sample_info = None
    for sample in samples:
        if sample['sample_name'] == sample_name:
            sample_info = sample
            break

    assert sample_info != None, "Something's wrong, I don't think that sample exists (%s)" % sample_name

    sample_info['library_name'] = config['general_params']['library_name']
    sample_info['sequencing_center'] = config['general_params']['sequencing_center']

    return read_group_string_template % sample_info
