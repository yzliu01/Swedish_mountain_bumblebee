#!/bin/bash -l

#SBATCH -A snic2021-22-195
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 08-00:00:00
#SBATCH --array=1-209%10
#SBATCH -J B.syl.all.sorted_bwa_mapping_SE
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu@imbim.uu.se
#SBATCH --output=bwa_mapping_all_sorted_SE_B.sylvicola.%A.o
#SBATCH --error=bwa_mapping_all_sorted_SE_B.sylvicola.%A.e

# load modules
module load bioinfo-tools bwa samtools bamtools

# directories to fastq files
SEQDIR=/x/x

# Reference sequence
REFDIR=/crex/proj/snic2020-6-58/private/yzliu/reference_genome
REF=$REFDIR/Bombus_sylvicola_v1.fna

# defining the SEED files
cd $SEQDIR
SEEDFILE1=seed1_R1_new				### Read1.fastq.gz
SEEDFILE2=seed2_R2_new				### Read2.fastq.gz
SEEDFILE3=seed3_RG					### ReadGroup for bwa
SEEDFILE4=seed4_sampleName_new		### Sample name

SEED1=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SEEDFILE1)
SEED2=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SEEDFILE2)
SEED3=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SEEDFILE3)
SEED4=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SEEDFILE4)

File1=${SEED4/.fastq.gz/.bam}	     		 ### bam file
File2=${File1/.bam/_sort.bam}			     ### sorted bam file

# application run with 8 threads

echo "run mapping for $SEED1 | $SEED2"

bwa mem -t 8 -R $SEED3 $REF $SEED1 $SEED2 | samtools view -b -@ 8 -o $File1 

samtools sort -@ 8 -m 5G -o $File2 $File1

samtools index $File2

bamtools stats -in $File2

