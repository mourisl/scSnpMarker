import sys
import argparse
import pysam

parser = argparse.ArgumentParser(description="Find SNPs that are unique to group 1's BAM files/barcode")
parser.add_argument("-l", help="list of BAM files with group ID", dest="list", required=True)
parser.add_argument("-g", help="list of genes and intervals in bed format (chr start end gene_id)", dest="gene", required=True)

args = parser.parse_args()

# 
geneIntervals = {}
fp = open(args.gene, "r")
for line in fp:
  cols = line.rstrip().split()
  gene = cols[3]
  if gene not in geneIntervals:
    geneIntervals[gene] = []
  geneIntervals[gene].append([cols[0], int(cols[1]), int(cols[2])]) # Coordinate is bed format
fp.close()

# collect the samples
bams = {} # key: bam file name, value: group
fp = open(args.list)
for line in fp:
  cols = line.rstrip().split()
  bams[cols[0]] = int(cols[1])
fp.close()

# Process each intervals
for g in geneIntervals:
  for label in [0, 1]:
    for b in bams:
      print(g, b)  
