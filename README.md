## Created bed file from gencode gtf
awk '$3=="exon"' gencode.v35.annotation.gtf | awk '{print $1"\t"$4-4"\t"$5"\t"$16}'| sed -e "s/[\";]//g" | sort -k1,1 -k2,2n > gencode_v35.bed
bedtools merge -c 4 -o distinct -i gencode_v35.bed > gencode_v35_merged.bed
