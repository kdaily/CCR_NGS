"""Dummy helper functions.

"""

import os

from ccrngspy import utils

logger = utils.make_local_logger("Dummy helper logging", level="debug", color="green")

def make_dummy_param_list(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    Sets the input and output to be the same.
    
    """

    final_list = []

    fastq_dir = config['general_params']['fastq_input_dir']
    output_dir = config['fastqc_params']['output_dir']
    
    for sample in samples:

        params = dict(sample=sample['samplename'])

        tmp1 = [os.path.join(fastq_dir, sample['filename1']),
                os.path.join(fastq_dir, sample['filename1']),
                params]
        tmp2 = [os.path.join(fastq_dir, sample['filename2']),
                os.path.join(fastq_dir, sample['filename2']),
                params]
        
        final_list.extend([tmp1, tmp2])

    return final_list
