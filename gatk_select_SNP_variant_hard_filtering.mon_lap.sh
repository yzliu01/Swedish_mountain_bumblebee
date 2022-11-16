#!/bin/bash -l

#SBATCH -A snic2021-22-195
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:30:00
##SBATCH --array=1-592%100
#SBATCH -J select_SNP_hard_filtering_89haploid_118diploid_B_syl
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu@imbim.uu.se
#SBATCH --error=select_SNP_hard_filtering_89haploid_118diploid_B_syl.%A.e
#SBATCH --output=select_SNP_hard_filtering_89haploid_118diploid_B_syl.%A.o

# Main directories as variables
# Path to where you have the genotyped vcf-files
SEQDIR=/x/x

# Reference sequence
# Path and file name for ref genome
REF=/x/reference_genome/Bombus_sylvicola_v1.fna

cd $SEQDIR

INTPUT_NAMES=genotyping_output_89haploid_118diploid_B_syl.g.vcf.gz
#OUTPUT_NAMES=${INTPUT_NAMES/.g.vcf.gz/.SNP.vcf.gz} # modify file names

module load bioinfo-tools
module load GATK

# select SNPs
gatk SelectVariants -R $REF \
-V $INTPUT_NAMES \
--select-type-to-include SNP \
-O gatk_89haploid_118diploid_B_syl.SNP.vcf.gz

# hard filtering
gatk VariantFiltration \
-R $REF \
-V gatk_89haploid_118diploid_B_syl.SNP.vcf.gz \
-O gatk_89haploid_118diploid_B_syl.SNP_hard_filtered.vcf.gz \
--filter-name "QD" \
--filter-expression "QD<2.0" \
--filter-name "FS" \
--filter-expression "FS>60.0" \
--filter-name "MQ" \
--filter-expression "MQ<40.0" \
--filter-name "MQRankSum" \
--filter-expression "MQRankSum<-12.5" \
--filter-name "ReadPosRankSum" \
--filter-expression "ReadPosRankSum<-8.0" \
--filter-name "SOR" \
--filter-expression "SOR>3.0"
