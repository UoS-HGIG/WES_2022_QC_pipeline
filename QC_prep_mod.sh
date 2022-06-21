#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#run this from the directory that your joint-called vcf is in
cd $SLURM_SUBMIT_DIR

## $1= tab delimited mapping file of of samples and their matching directories for a batch e.g.:
# AD001	/mainfs/hgig/private/EXOME_DATA/Macrogen_feb_2021/AD001
# AD023	/mainfs/hgig/private/EXOME_DATA/Macrogen_feb_2021/AD023
# AD005	/mainfs/hgig/private/EXOME_DATA/Macrogen_feb_2021/AD005
# AD020	/mainfs/hgig/private/EXOME_DATA/Macrogen_feb_2021/AD020
# AD018	/mainfs/hgig/private/EXOME_DATA/Macrogen_feb_2021/AD018


while read -r line; do
   ID=$(echo "$line"| cut -f1)
   directory=$(echo "$line"| cut -f2)
   cd $directory
   cp /mainfs/hgig/private/EXOME_DATA/PIPELINE_SCRIPTS_I5/QC/runbedtools.sh .
   sbatch runbedtools.sh $ID
   cd $SLURM_SUBMIT_DIR/
done < $1

#Will result in output of ${sample}_coverage.stats and ${sampleID}_coverage.hist` in individual sample directories

cp /mainfs/hgig/private/EXOME_DATA/PIPELINE_SCRIPTS_I5/QC/automateQC_shiny_modified.sh .

#When all jobs complete, run automateQC.sh from directory with joint called sample in
# sh automateQC_shiny_modified.sh <joint_called_dataset.vcf.gz> <output> <mapping_file>
