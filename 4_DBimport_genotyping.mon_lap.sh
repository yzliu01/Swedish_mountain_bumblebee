#!/bin/bash -l

#SBATCH -A snic2021-22-195
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 9-00:00:00
#SBATCH -a 1-592
#SBATCH -J B.syl_DBimport_genotyping
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu@imbim.uu.se
#SBATCH --error=B.syl_genotyping.%A.e
#SBATCH --output=B.syl_genotyping.%A.o

# path to vcf files
SEQDIR=/x/x

# Reference sequence
REFDIR=/crex/proj/snic2020-6-58/private/yzliu/reference_genome
REF=$REFDIR/Bombus_sylvicola_v1.fna

# The sample-name-map that you created in the previous step (tab-separated list of each sample name followed by the full path to its .g.vcf-file. One line per sample.)
# example:
# M_L001	/crex/proj/snic2020-6-58/private/seq_data/P23261_feb22_bombus_osmia/PAS/gVCF/M_L001.sorted_marked_dups.bam.g.vcf.gz
SAMP_NAME=$SEQDIR/sample_name_map_B.sylvicola.new.txt

        # The path to the directory where you want to have the output files
        # I'd recommend you create a separate directory for this

        # Create directory for DataBase and VCF fiels
		cd $SEQDIR; mkdir DB_VCF; cd DB_VCF; mkdir vcfs_all_indv DB_scaffold
		
OUTDIR=$SEQDIR/DB_VCF/vcfs_all_indv
DBDIR=$SEQDIR/DB_VCF/DB_scaffold

        # Get all the contig (or scaffold) names in array from the reference genome fasta file
scaffold_names=$(grep '>' $REF | sed 's/>//g' | sort -V | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Haploid or diploid genotyping?
# PLOIDY=2
# -ploidy $PLOIDY \

        module load bioinfo-tools
        module load GATK
        # Create a database for each contig/scaffold
        # The following command is submitted in parallel on multiple cores
echo "Running GenomicsDBImport for $scaffold_names"
gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
--sample-name-map $SAMP_NAME \
--genomicsdb-workspace-path $DBDIR/${scaffold_names} \
--intervals ${scaffold_names}

        # Use the databases for joint genotyping
        # The following command is submitted in array on multiple cores
echo "Running GenotypeGVCFs for ${scaffold_names}"
gatk --java-options "-Xmx4g -Xms4g" GenotypeGVCFs \
-R $REF \
-V gendb://$DBDIR/${scaffold_names} \
-new-qual true \
-O $OUTDIR/${scaffold_names}_output.g.vcf.gz


