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


