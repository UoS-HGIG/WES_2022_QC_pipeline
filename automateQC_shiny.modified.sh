#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --nodes=1
cd $SLURM_SUBMIT_DIR

module load picard
module load GATK/4.1.2
module load jdk/8.u181
module load biobuilds
module load R

DBSITES=/mainfs/hgig/public/HUMAN_REFS/HG38/common_dbsnp_151.hg38.vcf
REF=/mainfs/hgig/public/HUMAN_REFS/HG38/REF_HLA/GRCh38_full_analysis_set_plus_decoy_hla.fa

##IMOGEN MODIFY HERE##
CAPTURE=/mainfs/hgig/public/HUMAN_REFS/HG38/target_capture/SSV5_6_union_gc_pad150bp.bed

#First argument $1 is input
#Second argument $2 is output
#third argument $3 is mapping file (tab delimited file of sample ID and directory path)

##gatk beta version
 gatk VariantEval \
        -R ${REF} \
           -O ${2}.QC \
           --eval $1

##picard
picard CollectVariantCallingMetrics DBSNP=${DBSITES} \
    I=$1 \
    O=$2

#output from GATK is: SAMPLE.QC
#output from picard is: SAMPLE.variant_calling_detail_metrics and SAMPLE.variant_calling_summary_metrics

#deal with GATK_output
grep -v '#' ${2}.QC > tmp
grep 'CompOverlap' tmp > CompOverlap.sp
grep 'CountVariants' tmp > CountVariants.sp
grep 'IndelSummary' tmp > IndelSummary.sp
grep 'MetricsCollection' tmp > MetricsCollection.sp
grep 'MultiallelicSummary' tmp > MultiallelicSummary.sp
grep 'TiTvVariantEvaluator' tmp > TiTvVariantEvaluator.sp

find $PWD -name \*.sp > tsv.list

#replace spaces with tabs
while read -r file; do sed 's/ \+ /\t/g' $file > ${file}.tsv; done < tsv.list

find $PWD -name \*.sp.tsv > tsv.list

#remove unwanted columns and rows
while read -r file; do cut -f1,6- $file > ${file}_updated; head -2 ${file}_updated > ${file}_table; done < tsv.list

#concatenate columns
paste MetricsCollection.sp.tsv_table CompOverlap.sp.tsv_table CountVariants.sp.tsv_table TiTvVariantEvaluator.sp.tsv_table MultiallelicSummary.sp.tsv_table IndelSummary.sp.tsv_table > QC.${2}.output

mkdir QC_output_for_shiny
cp MetricsCollection.sp.tsv_table CompOverlap.sp.tsv_table CountVariants.sp.tsv_table TiTvVariantEvaluator.sp.tsv_table MultiallelicSummary.sp.tsv_table IndelSummary.sp.tsv_table ./QC_output_for_shiny

#deal with histogram separately
grep 'IndelLengthHistogram' tmp > IndelLengthHistogram.sp
sed 's/ \+ /\t/g' IndelLengthHistogram.sp > IndelLengthHistogram.sp_updated
cut -f1,6- IndelLengthHistogram.sp_updated | head -21 > ./QC_output_for_shiny/histogram.tsv

#deal with picard detail_metrics
for file in *detail_metrics; do grep -v '#' $file > ${file}.table; done

#deal with picard summary metrics
for file in *summary_metrics; do grep -v '#' $file > ${file}.table; done

cp ${2}.variant_calling_detail_metrics.table ./QC_output_for_shiny/joint_called.variant_calling_detail_metrics.table
cp ${2}.variant_calling_summary_metrics.table ./QC_output_for_shiny/joint_called.variant_calling_summary_metrics.table
cp QC.${2}.output ./QC_output_for_shiny/joint_called.QC

rm *_updated *.list tmp *.sp *.sp.tsv

#Make sure runbedtools.sh has finished and outputted coverage stats and histograms
#find $PWD -name \*_coverage.hist > coverage.hist.paths
#find $PWD -name \*_coverage.stats > coverage.stats.paths

cut -f2 mapping_file_shiny > b
while read -r line; do ls ${line}*coverage*hist; done < b > coverage.hist.paths
while read -r line; do ls ${line}*coverage*stats; done < b > coverage.stats.paths

while read -r line; do cp $line ./QC_output_for_shiny ; done < coverage.hist.paths
while read -r line; do cp $line ./QC_output_for_shiny ; done < coverage.stats.paths

#rm coverage.hist.paths coverage.stats.paths

cp /mainfs/hgig/private/EXOME_DATA/PIPELINE_SCRIPTS_I5/QC/app.R ./QC_output_for_shiny
cp /mainfs/hgig/private/EXOME_DATA/PIPELINE_SCRIPTS_I5/QC/report* ./QC_output_for_shiny 

#run peddy in local environment
#peddy -p 4 --sites hg38 --plot --prefix your_sample AD_jc_targ_vqsr.vcf.gz AD.ped

###################################################################
# MAKE SURE YOU PUT PEDDY FILES IN QC_output_for_shiny directory! #
###################################################################

#file ready to run app.R (run on local machine)
