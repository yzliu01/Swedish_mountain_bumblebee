#!/bin/bash -l

#SBATCH -A snic2021-22-195
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 20:30:00
#SBATCH --array=1-209%105
#SBATCH -J mark_dup_1-209_index_B.syvicola
#SBATCH --mail-type=all
#SBATCH --error=picard_samtools_mark_dups_index_1-209_B.syvicola.%A.e
#SBATCH --output=picard_samtools_mark_dups_index_1-209_B.syvicola.%A.o
#SBATCH --mail-user=yuanzhen.liu@imbim.uu.se

# directory to bam-files
SEQDIR=/x/x

# load modules
module load bioinfo-tools picard samtools

cd $SEQDIR
INPUT_NAMES=$(ls *_sort.bam | sed -n ${SLURM_ARRAY_TASK_ID}p) ## e.g. P23261_469_sort.bam
LINE_NUMBER=${SLURM_ARRAY_TASK_ID}
#OUTPUT_NAMES=$(awk '{if(NR==$LINE_NUMBER) print $5}' seed3_ReadGroup) ## M_L001
OUTPUT_NAMES=$(awk -v var=$LINE_NUMBER '{if(NR==var) print $5}' seed3_ReadGroup)

echo -e ${INPUT_NAMES[*]}"\t"${OUTPUT_NAMES[*]} >> SampleID_2_Name.1-209.log

java -Xmx8g -jar "/sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar" \
MarkDuplicates \
INPUT=$INPUT_NAMES \
OUTPUT=$OUTPUT_NAMES.sorted_marked_dups.bam \
METRICS_FILE=$OUTPUT_NAMES.sorted_marked_dups.csv >& $OUTPUT_NAMES.sorted_marked_dups.log
# output file: M_L001.sorted_marked_dups.bam

# Index the file
samtools index $OUTPUT_NAMES.sorted_marked_dups.bam
samtools idxstats $OUTPUT_NAMES.sorted_marked_dups.bam > $OUTPUT_NAMES.sorted_marked_dups.idxstats.csv

