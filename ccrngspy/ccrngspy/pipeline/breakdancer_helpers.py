import csv
import os
import sys

class BreakdancerParser:
    """Parse text output of Breakdancer.

    """

    def __init__(self, fh):
        self.fh = fh
        self.software = self.fh.readline()
        self.command = self.fh.readline()
        tmp = self.fh.readline()
        self.library_stats = self.fh.readline()

        self.reader = csv.DictReader(self.fh, delimiter="\t")

    def __iter__(self):
        return self

    def next(self):
        return self.reader.next()

def filter_by_chromosome(x, chrom):
    return filter(lambda y: y['Chr1'] == chrom or y['Chr2'] == chrom, x)

def process(directory, chrom, oh=sys.stdout):
    flist = filter(lambda x: x.endswith("txt"), os.listdir(directory))

    results = map(lambda x: (x, filter(lambda y: y['Chr1'] == chrom or y['Chr2'] == chrom, BreakdancerParser(file(x)))), flist)

    fieldnames = ['filename', 'Chr1', 'Pos1', 'Chr2', 'Pos2', 'num_Reads', 'Score', 'Size']
    wr = csv.DictWriter(oh, fieldnames=fieldnames, delimiter="\t")

    wr.writerow(dict(zip(fieldnames, fieldnames)))

    for (fname, rec) in results:
        for result in rec:
            row = dict(filename=fname, Chr1=result['Chr1'], Pos1=result['Pos1'], 
                       Chr2=result['Chr2'],  Pos2=result['Pos2'],
                       num_Reads=result['num_Reads'], 
                       Score=result['Score'], Size=result['Size'])

            wr.writerow(row)

def process_bed(directory, chrom, oh=sys.stdout):
    flist = filter(lambda x: x.endswith("txt"), os.listdir(directory))

    results = map(lambda x: (x, filter(lambda y: y['Chr1'] == chrom or y['Chr2'] == chrom, BreakdancerParser(file(x)))), flist)

    fieldnames = ['chrom', 'start', 'stop', 'name', 'score', 'strand']
    wr = csv.DictWriter(oh, fieldnames=fieldnames, delimiter="\t")

    for (fname, rec) in results:

        if rec:
            oh.write("track name='%s' description='%s'\n" %(fname, fname))

            for result in rec:
                row = dict(chrom=result['Chr1'], start=result['Pos1'],  stop=int(result['Pos1']) + 1,
                           name=".", score="0", strand="+")
                wr.writerow(row)
