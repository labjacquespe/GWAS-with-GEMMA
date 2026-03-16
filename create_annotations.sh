 #!/bin/sh

ml bedtools 

#----------------------------------
### 1) Download annotations file

# Gene annotations:
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz

# SNP annotations:
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
mkdir annotations
mv Arabidopsis_thaliana.TAIR10.62.gff3.gz annotations/.
mv 1001genomes_snp-short-indel_only_ACGTN.vcf.gz annotations/.


#----------------------------------
### 2) Extract SNPs and remove INDEL
zcat annotations/1001genomes_snp-short-indel_only_ACGTN.vcf.gz | \
 awk 'length($4)==1 && length($5)==1 {print $1":"$2"\t"$4"\t"$5}' \
 sort -k1b,1 > annotations/1001genomes_snp-short-indel_only_ACGTN.sorted.tsv


#----------------------------------
### 3) Extract SNP included in dataset and overlap with gene regions and intergenic regions.
zcat data/aragwas.genotypes_384.homozygote.vcf.gz | grep -v "#" | awk '{OFS="\t"; print $1,$2-1,$2,$3}' \
  > annotations/aragwas.genotypes_384.homozygote.bed
gzip annotations/aragwas.genotypes_384.homozygote.bed

#----------------------------------
### 4) Extract gene region annotations

### Extract gene regions
zcat annotations/Arabidopsis_thaliana.TAIR10.62.gff3.gz | awk '$3=="gene" && $1!="Mt" && $1!="Pt" {OFS="\t"; split($9,n,";"); print $1,$4,$5,"gene",n[1],n[2],n[3],n[4],$7}' |\
  sed 's/ID=gene://g' | sed 's/Name=//g' |sed 's/biotype=//g' | sed 's/description=//g' | \
  awk '{OFS="\t"; if($9=="+"){print $1,$2,$3,$4,$5,$2,$6,$7,$8,$9} if($9=="-") {print $1,$2,$3,$4,$5,$3,$6,$7,$8,$9}}' \
  > annotations/Arabidopsis_thaliana.TAIR10.62.gene.bed

### Intersect with SNPs and calculate TSS distance
bedtools intersect -wo \
  -a annotations/aragwas.genotypes_384.homozygote.bed.gz \
  -b annotations/Arabidopsis_thaliana.TAIR10.62.gene.bed | \
  awk '{OFS="\t"; dist=$10-$3; debut=$6-$3; fin=$7-$3; len=$7-$6; if (dist <0) {print $4,$8,$9,dist,len,debut,fin,$11,$12,$13,$14,"NA","NA","NA","NA","NA","NA","NA","NA","NA";next} {print $4,$8,"NA","NA","NA","NA","NA","NA","NA","NA","NA",$9,dist,len,debut,fin,$11,$12,$13,$14}}' | \
  sort -k1b,1 > annotations/overlap.gene.bed

### Extract gene component annotations
zcat annotations/Arabidopsis_thaliana.TAIR10.62.gff3.gz | \
  awk '$2=="araport11" && $3!="gene" && $3!="mRNA" && $3!="exon" && $1!="Mt" && $1!="Pt" {OFS="\t"; print $1,$4,$5,$3}' \
  > annotations/Arabidopsis_thaliana.TAIR10.62.gene_annotation.bed

### Intersect annotations with SNPs
bedtools intersect -wo \
  -a annotations/aragwas.genotypes_384.homozygote.bed.gz \
  -b annotations/Arabidopsis_thaliana.TAIR10.62.gene_annotation.bed | \
  cut -f 4,8 | sort -k1b,1 |\
  uniq > annotations/overlap.gene_annotation.bed

#----------------------------------
### 5) Intergenic regions

### Complement of gene regions
bedtools complement -i annotations/Arabidopsis_thaliana.TAIR10.62.gene.bed \
  -g annotations/genome.txt | \
  bedtools closest -a - -b annotations/Arabidopsis_thaliana.TAIR10.62.gene.bed -D ref | \
  awk '{OFS="\t"; print $1,$2,$3,"intergenic",$8,$9,$5,$6,$10,$11,$12,$13}' \
  > annotations/Arabidopsis_thaliana.TAIR10.62.intergenic.bed

### Intersect with SNPs and calculate TSS distance
bedtools intersect -wo \
 -a annotations/aragwas.genotypes_384.homozygote.bed.gz \
 -b annotations/Arabidopsis_thaliana.TAIR10.62.intergenic.bed | \
 awk '{OFS="\t"; dist=$10-$3;len=$12-$11;debut=$11-$3;fin=$12-$3; print $4,$8,$9,dist,len,debut,fin,$13,$14,$15,$16}' | awk '{key = $1; $1 = "" ; sub(/^ /, "", $0);
  if (key in a)
    a[key] = a[key] "\t" $0
  else
    a[key] = $0} END {for (k in a) {print k "\t" a[k]} }'  | tr ' ' '\t' | sort -Vk1,1 > data/annotations/overlap.intergenic.concat.bed

### end of chromosomes (n=1,378)
awk 'NF==11 {OFS="\t"; if ($4<=0) {print $0,"NA","NA","NA","NA","NA","NA","NA","NA","NA"} if ($4>0) {print $1,$2,"NA","NA","NA","NA","NA","NA","NA","NA","NA",$3,$4,$5,$6,$7,$8,$9,$10,$11}}' \
  annotations/overlap.intergenic.concat.bed \
  > annotations/overlap.intergenic.harmonized.bed

### one each side (n=4,683,163)
awk 'NF==21' annotations/overlap.intergenic.concat.bed | cut -f 1-11,13- \
  >> annotations/overlap.intergenic.harmonized.bed

### 2 genes with same start (n=2,548)
awk 'NF==31' annotations/overlap.intergenic.concat.bed | awk '{OFS="\t"; if($14<0) {if($14>$4) {print $1,$2,$13,$14,$15,$16,$17,$18,$19,$20,$21,$23,$24,$25,$26,$27,$28,$29,$30,$31; next} if($14<$4) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$23,$24,$25,$26,$27,$28,$29,$30,$31; next}} if($14>0) {if($14<$24) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$13,$14,$15,$16,$17,$18,$19,$20,$21; next} if ($11>$18) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$23,$24,$25,$26,$27,$28,$29,$30,$31}} if ($14==$24) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$13","$23,$14,$15","$25,$16","$26,$17","$27,$18","$28,$19","$29,$20","$30,$21","$31}}' \
  >> data/annotations/overlap.intergenic.harmonized.bed

#----------------------------------
### 6) Bring everything together
join -a1 annotations/overlap.gene.bed annotations/overlap.gene_annotation.bed | \
  tr ' ' '\t' | \
  cat - annotations/overlap.intergenic.harmonized.bed | \
  awk ' {if(NF==20) {print $0"\tNA"} if(NF==21) {print}}' | sort -k1b,1 | \
  join -a1 - annotations/1001genomes_snp-short-indel_only_ACGTN.sorted.tsv | \
  tr ' ' '\t' | uniq | \
  awk ' {if(NF==21) {print $0"\tNA\tNA"} if(NF==23) {print}}' | \
  gzip -c > annotations/annotations.len.tsv.gz



