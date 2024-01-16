#import sys
import argparse
import pysam

# Parse the ">>,,...ACGT"-like mpileup sequence  
def ParsePileupSequence(seq):
  nucCnt = {"A":0, "C":0, "G":0, "T":0, ".":0} # "." is for ref match
  for c in seq:
    c = c.upper()
    if (c not in "ACGT.," or c == ""):
      continue
    if (c == '.' or c == ','):
      nucCnt['.'] += 1
    else:
      nucCnt[c] += 1
  return [nucCnt[c] for c in "ACGT."]

parser = argparse.ArgumentParser(description="Find SNPs that are unique to group 1's BAM files/barcode")
parser.add_argument("-l", help="list of BAM files with group ID", dest="list", required=True)
parser.add_argument("-g", help="list of genes and intervals in bed format (chr start end gene_id)", dest="gene", required=True)
parser.add_argument("--min-cov", help="the minimum number of reads covering a base for a valid cell", dest="minCov", default=3)
parser.add_argument("--min-cell", help="the minimum number of cells that has valid base calling results in a group", dest="minCellSupport", default=3)
parser.add_argument("--debug", help="output the count information at each position", action="store_true")

args = parser.parse_args()

mismatchFrac = 0.99
minCov = args.minCov
minCellSupport = args.minCellSupport
groupCellFrac = [1, 0.6]

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

# Create the pysam object for each bam
pysams = {}
for b in bams:
  pysams[b] = pysam.AlignmentFile(b, "rb") 

# Process each intervals
for g in geneIntervals:
  for interval in geneIntervals[g]:
    cellSupportResult = [] 
    intervalLen = interval[2] - interval[1]
    for i in range(intervalLen):
      cellSupportResult.append([[0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]) # the first index is for label 0,1. the second index is for alternative nucleotide (the value [0..3], and the last is for )
    
    # Collect cell support
    for label in [0, 1]: # group label
      for b,bl in bams.items():
        if (bl != label):
          continue
        for pileupread in pysams[b].pileup(interval[0], interval[1], interval[2] - 1, min_mapping_quality=1):
          if (pileupread.reference_pos < interval[1] or pileupread.reference_pos >= interval[2]):
            continue
          mpileupseq = pileupread.get_query_sequences()
          nucResult = ParsePileupSequence(mpileupseq) 
          totalResultCount = sum(nucResult)

          #if (args.debug is not None):
          #  print(b, bl, pileupread.reference_pos, nucResult)

          if (totalResultCount < minCov):
            continue
          for j in range(5):
            if (nucResult[j] >= totalResultCount * mismatchFrac):
              cellSupportResult[pileupread.reference_pos - interval[1]][label][j] += 1
    
    # Go through each position to output the group 1 specific SNPs
    for i in range(intervalLen):
      support = cellSupportResult[i]
      if (args.debug):
        print(gene, interval[0], interval[1] + i, support)
      if (sum(support[0]) < minCellSupport or sum(support[1]) < minCellSupport):
        continue
      nucType = [-1, -1]
      for j in [0, 1]:
        totalCell = sum(support[j])
        for t, cnt in enumerate(support[j]):
          if (cnt > totalCell * groupCellFrac[j] - 1e-6): 
            nucType[j] = t
      if (nucType[0] != -1 and nucType[1] != -1 and nucType[0] != nucType[1]):
        print("%s %s %d %d %d %d %d"%(gene, interval[0], interval[1] + i, nucType[0], nucType[1], support[0][nucType[0]], support[1][nucType[1]]))
      
