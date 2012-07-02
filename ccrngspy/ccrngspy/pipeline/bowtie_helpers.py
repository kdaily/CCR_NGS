"""Helper files for pipeline bowtie steps.

"""

import os

from ccrngspy import utils

logger = utils.make_local_logger("Bowtie helper logging", level="debug", color="green")

def make_bowtie_param_list(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    Needs to be a list of [[input1, input2], output, params]; for the fastqc script.

    The output is the file RUM.sam.

    while the params are taken from the global opts variable
    (and possibly from the YAML config file).
    
    """

    final_list = []

    fastq_dir = config['general_params']['fastq_input_dir']
    log_dir = config['general_params']['log_file_dir']
    rum_dir = config['bowtie_params']['output_dir']

    
    for sample in samples:
        params = dict(sample=sample['samplename'])
    
        tmp = [[os.path.join(fastq_dir, sample['filename1']),
                os.path.join(fastq_dir, sample['filename2'])],
               os.path.join(bowtie_dir, sample['samplename'], "bowtie2.sam"),
               params]
        
        final_list.append(tmp)

    return final_list
