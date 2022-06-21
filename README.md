# WES_2022_QC_pipeline

############################################################################
Step1: run QC_prep_mod.sh mapping_file_shiny in directory with joint-called vcf
############################################################################

NB: first argument is tab delimited mapping file of sample ID and directory, e.g.

SAMPLE1 /path/to/SAMPLE1
SAMPLE2 /path/to/SAMPLE2

When all jobs complete (that is to say ${sample}_coverage.stats and ${sampleID}_coverage.hist files are generated for each sample in batch)... move to step 2

###########################################################################################
Step2: run automateQC_shiny_modified.sh <input_joint_called.vcf.gz> <output> <mapping_file>
###########################################################################################

NB: first argument = joint-called vcf of samples, second argument = desired output name, third argument = mapping file you used in QC_prep_mod.sh script
You run above from the directory that contains your joint-called vcf (the automateQC_shiny_modified.sh file will already be in here as it is copied across in QC_prep.sh)

###############################################################################################################
Step3: run peddy in your conda environment and save output to the newly created directory /QC_output_for_shiny
###############################################################################################################

NB you need a .ped file with the following headers:
Family ID       Individual ID   Paternal ID     Maternal ID     Sex (1=male; 2=female; other=unknown)   Phenotype

You run peddy using the following command:
peddy -p 4 --sites hg38 --plot --prefix SAMPLES joint_called_targeted+150bp.recalibrated_VQSR.vcf.gz joint_called.ped

###############################################################
Step4: download QC_output_for_shiny directory to local machine
###############################################################

######################################################################################
Step5: open Rstudio and make the /QC_output_for_shiny directory your working directory
######################################################################################

Make sure you have all the libraries installed in app.R

########################################
Step6: in Rstudio console run `runApp()` 
########################################
