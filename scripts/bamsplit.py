#!/usr/bin/env python
# split a bam file by read group ID
# Sean Davis <seandavi@gmail.com>
# March 10, 2012
# From: https://gist.github.com/kdaily/4955971

import pysam
import argparse
import logging
import os

logging.basicConfig(level=logging.INFO)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('filename',help="The bam filename")
    parser.add_argument('--output_dir', help="Output directory", default=".")
    parser.add_argument('--output_file_expr', help="Format for new filenames", default=r"%(ID)s_%(PL)s_%(LB)s_%(PU)s_%(SM)s.bam")
 
    opts = parser.parse_args()
 
    infile = pysam.Samfile(opts.filename,'rb')
 
    header = infile.header

    readgroups = header['RG']

    # remove readgroups from header
    del(header['RG'])

    outfiles = {}

    for rg in readgroups:
        tmphead = header
        tmphead['RG'] = [rg]
        bam_out_filename = os.path.join(opts.output_dir, opts.output_file_expr % rg)

        logging.info('Creating new BAM file: %s', (bam_out_filename, 'wb'))
        outfiles[rg['ID']] = pysam.Samfile(bam_out_filename, 'wb', header=tmphead)

    j = 0
    for read in infile:
        j += 1
        idtag = [x[1] for x in read.tags if x[0] == 'RG'][0]
        outfiles[idtag].write(read)

        if((j % 1000000) == 0):
            logging.info('read and wrote %d records', (j))
        

    for f in outfiles.values():
        f.close()

    infile.close()

if __name__ == "__main__":
    main()
