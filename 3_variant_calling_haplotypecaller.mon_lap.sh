#!/bin/bash -l

#SBATCH -A snic2021-22-195
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 3-10:00:00
#SBATCH -J Bsyl_haplotypecaller
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu@imbim.uu.se
#SBATCH --error=variant_calling_haplotypecaller_Bsyl.%A.e
#SBATCH --output=variant_calling_haplotypecaller_Bsyl.%A.o

# path to bam-files
SEQDIR=/x/x
OUTDIR=$SEQDIR"/gVCF"
cd $SEQDIR
mkdir gVCF

# Update with the path to your ref genome
REFDIR=/crex/proj/snic2020-6-58/private/yzliu/reference_genome
# Update this with the name of the ref fasta file
REF=$REFDIR/Bombus_sylvicola_v1.fna

# Haploid or diploid genotyping
PLOIDY=2

# This function is to run the jobs in parallel
function pwait() {
  while [ $(jobs -p | wc -l) -ge $1 ]; do
    sleep $2
  done
}

function call_variants() {
        module load bioinfo-tools
        module load GATK #Need to use GATK4 for GenomicsDBImport

        cd $SEQDIR

        samp_list=`ls *sorted_marked_dups.bam`
		# example: M_L010.sorted.marked_dups.bam
        for sample in ${samp_list[*]}; do
#		for sample in $samp_list; do
                # The jobs are submitted in the background (with the "&" option below)
                # The pwait function prevents jobs from being submitted if the maximum number of jobs (40) are already running
                # Then it waits for 10 s and then checks again
				
#                pwait 40 10s

                echo "Running HaplotypeCaller for $sample"

                gatk HaplotypeCaller \
                -R $REF \
                -I $sample \
                -O $OUTDIR/$sample.g.vcf.gz \
                -ERC GVCF \
                -ploidy $PLOIDY &

        done
        wait # wait until all jobs have finished

}

call_variants
