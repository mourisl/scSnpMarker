## Step 0.1 Create bam_group file
A two-column list file, the first column is the bam file, the second is the indicator: -1, 0, 1. -1 is for unknown condition cell bam file, 0 can be for normal cells, 1 can be for cancer cells (0,1 is just for condition purpose).

## Step 0.2 Created bed file from gencode gtf
awk '$3=="exon"' gencode.v35.annotation.gtf | awk '{print $1"\t"$4-4"\t"$5"\t"$16}'| sed -e "s/[\";]//g" | sort -k1,1 -k2,2n > gencode_v35.bed
bedtools merge -c 4 -o distinct -i gencode_v35.bed > gencode_v35_merged.bed

## Generate informative SNPs
python3 FindMarkerSnps.py -l bam_group -g gencode_v35_merged.bed > informative_snp.list

## Classifiy unknown samples
python3 ClassifyFromMarkerSnps.py -s informative_snp.list -l bam_group > classify_result.out
