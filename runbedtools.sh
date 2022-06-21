#!/bin/bash
#SBATCH --ntasks=40
#SBATCH --time=01:30:00
#SBATCH --nodes=1
cd $SLURM_SUBMIT_DIR

#########################################################################################
# MAPPING STATISTICS
#########################################################################################

## load modules ##
module load biobuilds/2017.05
module load samtools/1.9
module load jdk

#load variables
# First argument is sample ID

## coverage calculations ##
samtools view -bq 20 -F 1796 "$1".GATK.recal.bam | bamToBed -i stdin > final.bed
wc -l final.bed
awk '{if ($2 != $3) {print $0}}' final.bed > temp.bed && mv temp.bed final.bed 
wc -l final.bed

#sort final.bed
sort -k1,1 -k2,2n final.bed > sorted.bed

## Changed around file input (as of bedtools v2.24.0)

##IMOGEN MODIFY CAPTURE KIT HERE##
coverageBed -hist -sorted -a /mainfs/hgig/public/HUMAN_REFS/capture_kits/hg38_V456_union_QC_sorted.bed -b sorted.bed > coverage_gencode_interval.bed
coverageBed -sorted -a /mainfs/hgig/public/HUMAN_REFS/capture_kits/hg38_V456_union_QC_sorted.bed -b sorted.bed > coverage_agilent_targets.bed
coverageBed -sorted -a /mainfs/hgig/private/EXOME_DATA/joint_call_2022/intervals/hg38_V456_union_QC_sorted+150bp.bed -b sorted.bed > coverage_agilent_targets_150.bed

## efficiency of capture ##
reads=`awk 'END {OFS = "\t"; print NR}' sorted.bed`

mapped_to_target_reads=`awk '{SUM += $4} END {OFS = "\t";print SUM}' coverage_agilent_targets.bed`
# uncomment following if Ckit V5
#mapped_to_target_reads=`awk '{SUM += $5} END {OFS = "\t";print SUM}' coverage_agilent_targets.bed`
percent1=`awk 'BEGIN{printf("%0.2f", ('$mapped_to_target_reads' / '$reads') * 100)}'`

mapped_to_target_reads_plus_150=`awk '{SUM += $4} END {OFS = "\t";print SUM}' coverage_agilent_targets_150.bed`
# uncomment following if Ckit V5
#mapped_to_target_reads_plus_150=`awk '{SUM += $5} END {OFS = "\t";print SUM}' coverage_agilent_targets_150.bed`
percent2=`awk 'BEGIN{printf("%0.2f", ('$mapped_to_target_reads_plus_150' / '$reads') * 100)}'`

## coverage ##
grep all coverage_gencode_interval.bed > "$1"_coverage.hist
meancov=`awk '{if ($2>=1) (SUM += $2*$5)} END {printf ("%0.2f", SUM)}' "$1"_coverage.hist`

## completeness of coverage ##
cov1xpc=`awk '{if ($2>=1) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' "$1"_coverage.hist`
cov5xpc=`awk '{if ($2>=5) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' "$1"_coverage.hist`
cov10xpc=`awk '{if ($2>=10) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' "$1"_coverage.hist`
cov20xpc=`awk '{if ($2>=20) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' "$1"_coverage.hist`
cov0x=`awk '{if ($2>=0) (SUM += $3)} END {print SUM}' "$1"_coverage.hist`
cov1x=`awk '{if ($2>=1) (SUM += $3)} END {print SUM}' "$1"_coverage.hist`
cov5x=`awk '{if ($2>=5) (SUM += $3)} END {print SUM}' "$1"_coverage.hist`
cov10x=`awk '{if ($2>=10) (SUM += $3)} END {print SUM}' "$1"_coverage.hist`
cov20x=`awk '{if ($2>=20) (SUM += $3)} END {print SUM}' "$1"_coverage.hist`

#add sample ID to coverage.hist - step required for shiny app
for f in *coverage.hist; do sed -i "s/$/\t$f/" $f; done

#report generation
printf "sample\t"$1"\n" > "$1"_coverage.stats
printf "total_reads\t"$reads"\n" >> "$1"_coverage.stats
printf "mapped_to_target_reads\t"$mapped_to_target_reads"\n" >> "$1"_coverage.stats
printf "percentage\t"$percent1"\n" >> "$1"_coverage.stats
printf "mapped_to_target_reads_plus_150bp\t"$mapped_to_target_reads_plus_150"\n" >> "$1"_coverage.stats
printf "percentage\t"$percent2"\n" >> "$1"_coverage.stats
printf "mean_coverage\t"$meancov"\n" >> "$1"_coverage.stats
printf "accessible_target_bases\t"$cov0x"\n" >> "$1"_coverage.stats
printf "accessible_target_bases_1x\t"$cov1x"\n" >> "$1"_coverage.stats
printf "percentage\t"$cov1xpc"\n" >> "$1"_coverage.stats
printf "accessible_target_bases_5x\t"$cov5x"\n" >> "$1"_coverage.stats
printf "percentage\t"$cov5xpc"\n" >> "$1"_coverage.stats
printf "accessible_target_bases_10x\t"$cov10x"\n" >> "$1"_coverage.stats
printf "percentage\t"$cov10xpc"\n" >> "$1"_coverage.stats
printf "target_bases_20x\t"$cov20x"\n" >> "$1"_coverage.stats
printf "percentage\t"$cov20xpc"\n" >> "$1"_coverage.stats

rm final.bed
