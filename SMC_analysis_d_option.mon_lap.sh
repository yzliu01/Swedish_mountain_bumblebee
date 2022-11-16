#!/bin/bash -l

#SBATCH -A snic2022-22-258
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 20:30:00
#SBATCH -J mon_10F_lap_11F_depth_X15_smc++
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanzhen.liu@imbim.uu.se
#SBATCH --error=mon_10F_lap_11F_depth_X15_smc_disting.%A.e
#SBATCH --output=mon_10F_lap_11F_depth_X15_smc_disting.%A.o

module load bioinfo-tools vcftools SMC++/1.15.2 

# 2.Convert your VCF(s) to the SMC++ input format with vcf2smc
# composite likelihood to obtain good estimates
VCF=/x/x
OUT=/x/analysis/SMC
REF=/x/reference_genome
contig=`awk '{print $1}' $REF/Bombus_sylvicola_v1.fna.fai | sort -V`

cd $OUT

# This function is to run the jobs in parallel
	function pwait() {
	while [ $(jobs -p | wc -l) -ge $1 ]; do
	sleep $2
	done
	}

## use multiple variables in for loop
## https://stackoverflow.com/questions/11215088/bash-shell-script-two-variables-in-for-loop
function mon_smc_parallel() {

		pwait 100 10s
		
	# use multiple variables in for loop
	for i in $contig; do
		for j in {BH_16,AC-020,AC-021,AC-034,M_L210,M_L213,BH_07,BH_08,BH_09,M_L119}; do
		smc++ vcf2smc -d $j $j $VCF/mon_10F_lap_11F_depth_X15.headline.vcf.gz \
		$OUT/mon/"Bombus_monticola.$i.$j.txt" $i \
		Bombus_monticola:BH_16,AC-020,AC-021,AC-034,M_L210,M_L213,BH_07,BH_08,BH_09,M_L119
		done

	done
	wait
}
#mon_smc_parallel

function lap_smc_parallel() {

		pwait 100 10s
		
	# use multiple variables in for loop
	for i in $contig; do
		for j in {BH_11,M_L021,M_L138,M_L180,M_L194,AC-006,AC-138,WO_003,WO_006,WO_383,WO_448,WO_449,WO_703,M_L018}; do
		smc++ vcf2smc -d $j $j $VCF/mon_10F_lap_11F_depth_X15.headline.vcf.gz \
		$OUT/lap/"Bombus_lapponicus.$i.$j.txt" $i \
		Bombus_lapponicus:BH_11,M_L021,M_L138,M_L180,M_L194,AC-006,AC-138,WO_003,WO_006,WO_383,WO_448,WO_449,WO_703,M_L018
		done

	done
	wait
}
#lap_smc_parallel

## run two functions simutaneously
mon_smc_parallel & lap_smc_parallel

# 3.Fit the model using estimate
for f in {lap,mon}; do
		# run parallel jobs
		pwait 2 10s
		smc++ estimate -o $OUT/$f 3.60e-9 $OUT/$f/Bombus_*.txt
		done

# 4.Visualize the results using plot
smc++ plot $OUT/mon_lap.composite_likelihood.plot.pdf $OUT/mon/model.final.json $OUT/lap/model.final.json

#smc++ plot $OUT/pas_mon_plot.pdf $OUT/pas/model.pas.json $OUT/mon/mon_model.final.json

## rename files using for
#OUT=/crex/proj/snic2020-6-58/private/yzliu/swedish_bumblebee/analysis/SMC
#cd $OUT
#for i in {$OUT/mon/model.final.json,$OUT/lap/model.final.json}; do
#mv $i $i.old
#done 

