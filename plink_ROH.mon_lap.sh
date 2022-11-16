#!/bin/bash -l

#SBATCH -A snic2022-22-258
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:10:00
#SBATCH -J ROH_mon_lap
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu@imbim.uu.se
#SBATCH --error=ROH_mon_lap.%A.e
#SBATCH --output=ROH_mon_lap.%A.o

module load bioinfo-tools
module load vcftools
module load plink

# input dir
INPUT_DIR=/x/x
PLINK_DIR=$INPUT_DIR/plink_format
ROH_DIR=$INPUT_DIR/homozygosity

# vcf files
cd $INPUT_DIR
# mon
INPUT_FILE=mon.DP3_GQ20_biSNP_maf005_f_missing05.infoDP4_3134Excesshet_25.rm_non_var.mon_49F.vcf.gz
ped_FILE=mon.DP3_GQ20_biSNP_maf005_f_missing05.infoDP4_3134Excesshet_25.rm_non_var.mon_49F
bed_FILE=mon.DP3_GQ20_biSNP_maf005_f_missing05.infoDP4_3134Excesshet_25.rm_non_var.mon_49F
plink_ROH_FILE=plink_ROH_mon_49F

# lap
INPUT_FILE=lap.DP3_GQ20_biSNP_maf005_f_missing05.infoDP4_3134Excesshet_25.rm_non_var.lap_53F.vcf.gz
ped_FILE=lap.DP3_GQ20_biSNP_maf005_f_missing05.infoDP4_3134Excesshet_25.rm_non_var.lap_53F
bed_FILE=lap.DP3_GQ20_biSNP_maf005_f_missing05.infoDP4_3134Excesshet_25.rm_non_var.lap_53F
plink_ROH_FILE=plink_ROH_lap_53F

# fam file sample ID issue
vcftools --gzvcf $INPUT_FILE --plink --out $PLINK_DIR/$ped_FILE.vcftools
plink --file $PLINK_DIR/$ped_FILE --make-bed --recode 12 --out $PLINK_DIR/$bed_FILE.vcftools

# run again to make ped file
plink --vcf $INPUT_FILE --make-bed --allow-extra-chr --out $PLINK_DIR/$ped_FILE

### for result
plink --bfile $PLINK_DIR/$bed_FILE --out $ROH_DIR/$plink_ROH_FILE.10kb_50snp_no_het.group --homozyg --homozyg-snp 50 --homozyg-kb 10 \
--homozyg-density 50 --homozyg-window-het 0 --homozyg-window-missing 3 --allow-extra-chr --homozyg-group

plink --bfile $PLINK_DIR/$bed_FILE --out $ROH_DIR/$plink_ROH_FILE.10kb_50snp_no_het.group --homozyg --homozyg-snp 50 --homozyg-kb 10 \
--homozyg-density 50 --homozyg-window-het 0 --homozyg-window-missing 3 --allow-extra-chr --homozyg-group

