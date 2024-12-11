#import sys
import argparse
import pysam
import utils

parser = argparse.ArgumentParser(description="Based on the FindMarkerSnp output to classify unknown BAM files")

parser.add_argument("-l", help="list of BAM files with group ID", dest="list", required=True)
parser.add_argument("-s", help="list of marker SNPs", dest="snpList", required=True)
parser.add_argument("--min-cov", help="the minimum number of reads covering a base for a valid cell", dest="minCov", default=3)
parser.add_argument("-f", help="fraction of read support that needs to be categorize to group 1", dest="frac", default=0.99)
parser.add_argument("--min-mapq", help="the minimum MAPQ of an alignment", dest="minMapQ", default=30)

args = parser.parse_args() 
minCov = args.minCov
minMapQ = args.minMapQ
frac = args.frac 

def IsBlacklistGenes(g):
  prefixes = ["IGHV", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ", 
      "TRBV", "TRBJ", "TRAV", "TRAJ", "TRGV", "TRGJ", "TRDV", "TRDJ"]#, "HLA"]
  for p in prefixes:
    if (g.startswith(p)):
      return True
  return False

# collect the samples
bams = {} # key: bam file name, value: group
fp = open(args.list)
for line in fp:
  cols = line.rstrip().split()
  bams[cols[0]] = int(cols[1])
fp.close()

# Collec the marker snps
snpList = []
fp = open(args.snpList, "r")
for line in fp:
  cols = line.rstrip().split()
  for i in range(2, len(cols)):
    cols[i] = int(cols[i])
  snpList.append(cols)
fp.close()

# Go through each file and assign the cateogry
for b in bams:
  #if (bams[b] != -1):
  #  continue
  sam = pysam.AlignmentFile(b, "rb") 
  groupCount = [0, 0, 0, 0] # support group 0. suport group 1. neither support. no enough read

  group0snps = []
  group1snps = []
  for snp in snpList:
    zerocov = 1 # zero coverage, will will be skipped by pileup 
    if (IsBlacklistGenes(snp[0])):
      continue

    for pileupread in sam.pileup(snp[1], snp[2], snp[2] + 1, min_mapping_quality=minMapQ):
      if (pileupread.reference_pos != snp[2]):
        continue
      mpileupseq = pileupread.get_query_sequences()
      nucResult = utils.ParsePileupSequence(mpileupseq) 
      totalResultCount = sum(nucResult)
      zerocov = 0

      if (totalResultCount < minCov):
        groupCount[3] += 1
        continue
      
      if (nucResult[ snp[4] ] >= frac * totalResultCount - 1e-6
          or nucResult[snp[4] ] > 0.5 * totalResultCount):
        group1snps.append("_".join([str(s) for s in snp[0:4]]) + "_" + "/".join([str(nucResult[ snp[4] ]), str(totalResultCount)]))
        groupCount[1] += 1
      elif (nucResult[snp[3]] > 0.5 * totalResultCount):
        group0snps.append("_".join([str(s) for s in snp[0:4]]) + "_" + "/".join([str(nucResult[ snp[3] ]), str(totalResultCount)]))
        groupCount[0] += 1
      else:
        groupCount[2] += 1
    groupCount[3] += zerocov 
  groupCall = -1
  
  #if (groupCount[1] > 0):
  #  groupCall = 1
  #elif (groupCount[0] > 0):
  #  groupCall = 0
  if (groupCount[1] + groupCount[0] >= 4):
    if (groupCount[1] > groupCount[0]):
      groupCall = 1
    elif (groupCount[0] > groupCount[1]):
      groupCall = 0;
  
  countStr = "\t".join([str(x) for x in groupCount])
  print("%s\t%d\t%d\t%s\t%s\t%s\t"%(b, bams[b], groupCall, countStr, 
    ";".join(group0snps), ";".join(group1snps)))
