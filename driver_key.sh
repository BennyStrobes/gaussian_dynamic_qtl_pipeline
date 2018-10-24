#!/bin/bash
#SBATCH --time=20:00:00 --partition=broadwl --mem=14GB

###############################################################################
# Input Data
###############################################################################

# Directory created by "time_step_independent_qtl_pipelines" scripts
# Contains 1 file per sample with information on each test (variant, target region)
# Each file (sample) has the same number of lines (tests)
cht_input_file_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/time_step_independent_qtl_pipelines/wasp/cht_input_files/"

# File containing all of the target regions we are using. We are using this file to convert from gene positions to ensable id
target_region_input_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/time_step_independent_qtl_pipelines/wasp/target_regions/target_regions_cis_distance_50000_maf_cutoff_0.1_min_reads_100_min_as_reads_25_min_het_counts_5_merged.txt"


# cell line specific pcs
cell_line_specific_pc_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/covariates/cell_line_ignore_missing_principal_components_9.txt"

# Gene expression data for all samples
# Expression is RPKM transformed, then quantile normalized.
# Script will standardize each gene
total_expression_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/processed_total_expression/quantile_normalized_no_projection.txt"

# Dosage-based genotypes for all samples
genotype_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/genotype/YRI_genotype.vcf"

# Directory containing results for negative binomial dynamic qtl model
nb_dynamic_qtl_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/dynamic_qtl_pipelines/ipsc_data_te/qtl_results/"

# File_name of per-time-step eqtl results (t=0)
t_0_eqtl_results_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/time_step_independent_qtl_pipelines/wasp/cht_output/cht_results_cis_distance_50000_maf_cutoff_0.1_min_reads_100_min_as_reads_25_min_het_counts_5_num_pc_3_time_0_eqtl_results.txt"

# Directory containing chromHMM results
# Each cell line has its own file with suffix $cell_line_identifier'_15_coreMarks_mnemonics.bed.gz'
chrom_hmm_input_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/chrom_hmm/"

# Time step independent data
# Suffixes are:
#### 1. "$time_step"_eqtl_results.txt, containing pvalues for each time step, along with MAF, and dist-toTSS for our variant-gene pairs
#### 2. "$time_step"_efdr_thresh_.1_significant_egenes.txt giving list of significant variant gene pairs at this time step
time_step_independent_stem="/project2/gilad/bstrober/ipsc_differentiation_19_lines/time_step_independent_qtl_pipelines/wasp/cht_output/cht_results_cis_distance_50000_maf_cutoff_0.1_min_reads_100_min_as_reads_25_min_het_counts_5_num_pc_3_time_"

# Directory containing gsea data
gsea_data_dir="/project2/gilad/bstrober/tools/tools/gsea/data/"


# File containing conversions from ensamble ids to gene symbols
gencode_file="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/gencode.v19.annotation.gtf.gz"

# File containing gwas hits
gwas_catalog="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/gwas_catalog_v1.0.2-associations_e92_r2018-07-17.tsv"

# Directory containing gtex gwas hits
gtex_gwas_hits_dir="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/gtex_gwas_data/"


# File containing list of cardiomyopathy genes
cardiomyopathy_gene_list="/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/cardiomyopathy_disease_genes/cardiomyopathy_genes.txt"



###############################################################################
# Output directories (aasume all of these exist prior to starting analysis)
###############################################################################

# Root directory for this of all ipsc data based results
output_root="/project2/gilad/bstrober/ipsc_differentiation_19_lines/gaussian_dynamic_qtl_pipelines/"

# Directory containing necessary input files to qtl tests
input_data_dir=$output_root"input_data/"


# Directory containing text files with results from dynamic qtl analysis
qtl_results_dir=$output_root"qtl_results/"

# Directory containing visualization of results found in qtl_results_dir
qtl_visualization_dir=$output_root"qtl_visualization/"

# Output directory containing chromHMM enrichment results
chrom_hmm_enrichment_directory=$output_root"chromHMM_enrichment/"

# Output directory containing tissue specific chromHMM enrichment results
tissue_specific_chrom_hmm_enrichment_dir=$output_root"tissue_specific_chromHMM_enrichment/"

# Output directory containing visualizations of top dynanmic qtl visualizations
dynamic_qtl_visualization_dir=$output_root"dynamic_qtl_hit_visualization/"

# Output directory containing gene set enrichments
gene_set_enrichment_dir=$output_root"gene_set_enrichment/"

# Output directory containing gene set enrichments within Marios' gene sets
cardiomyopathy_gene_set_enrichment_dir=$output_root"cardiomyopathy_gene_set_enrichment/"

# Output directory containing gwas overlaps
gwas_overlap_dir=$output_root"gwas_overlap/"

# Output directory containing comparisons to time-step independent analysis
time_step_comparison_dir=$output_root"time_step_independent_comparison/"

# Output directory containing figures for manuscripts
manuscript_figures_dir=$output_root"manuscript_figures/"


###############################################################################
# Dynamic QTL Calling
###############################################################################


##########################################
# Step 1: Create joint_test_input_file
##########################################
# joint_test_input_file is a table with 1 row per sample
# each row has 3 columns:
#### 1. Sample id
#### 2. environmental variable
#### 3. Absolute directory to CHT input file for that sample
# NOTE: This script is very specific to our data
# Takes less than a minute to run
# $environmental_variable is a parameter that describes how we parameterize the environmental variable. So far, this is done with:
### 1. 'time_steps': raw time format
environmental_variable_form="time_steps"
joint_test_input_file=$input_data_dir"joint_test_input_file_"$environmental_variable_form".txt"
if false; then
python create_joint_test_input_file.py $cht_input_file_dir $joint_test_input_file $environmental_variable_form $pseudotime_predictions_3_file $pseudotime_predictions_4_file $pseudotime_predictions_5_file $cell_line_specific_pc_file
fi



##########################################
# Step 2: Prepare input files for dynamic qtl calling
##########################################
# How to hanlde genotypes in the model. Options currently include:
## 1. "round"
## 2. "dosage"
genotype_version="dosage"

# Stem for all output files
input_data_file=$input_data_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_input.txt"
if false; then
python prepare_gaussian_dynamic_qtl_input_files.py $joint_test_input_file $target_region_input_file $total_expression_file $genotype_file $genotype_version $input_data_file $environmental_variable_form
fi

##########################################
# Step 3: Run GLM/GLMM dynamic qtl modeling
##########################################

#############################################################
# Whether to run:
### 1. 'glm'
### 2. 'glmm'
model_version="glmm"
# Whether to include coviates or not
#### 1. 'none'
#### 2. 'pc1'
#### 3. 'pc1_2'
#### 4. 'pc1_3'
#### 5. 'pc1_4'
#### 6. 'pc1_5'
covariate_method="none"
num_jobs="10"

permute="False"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
    output_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_"$permute"_results_"$job_number".txt"
    sbatch run_gaussian_dynamic_qtl.sh $input_data_file $output_file $model_version $permute $covariate_method $job_number $num_jobs
done
fi



permute="True"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
    output_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_"$permute"_results_"$job_number".txt"
    sbatch run_gaussian_dynamic_qtl.sh $input_data_file $output_file $model_version $permute $covariate_method $job_number $num_jobs
done
fi






##########################################
# Step 3: Multiple testing correction and visualization
##########################################

model_version="glm_quadratic"
covariate_method="pc1_5"
real_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"
null_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_True_results.txt"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
if false; then
sh multiple_testing_correction_and_visualization.sh $real_file $null_file $qtl_visualization_dir $parameter_string $qtl_results_dir $num_jobs $dynamic_qtl_visualization_dir $input_data_file $model_version $covariate_method
fi



##########################################
# Step 4: Temporal QTL enrichment analysis
##########################################

covariate_method="pc1_5"

# All variant-gene pairs tested and their corresponding pvalues (using GLM)
gaussian_glm_qtl_results_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glm_covariate_method_"$covariate_method"_permute_False_results.txt"
# All significant variant gene pairs (eFDR <= .01) and their corresponding pvalues (using GLM)
gaussian_glm_qtl_significant_results=$qtl_results_dir$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glm_covariate_method_"$covariate_method"_efdr_.01_significant.txt"
# All significant genes (eFDR <= .01) and their top variant and their corresponding pvalues (using GLM)
gaussian_glm_qtl_significant_gene_results=$qtl_results_dir$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glm_covariate_method_"$covariate_method"_efdr_.01_significant_egenes.txt"
# All variant-gene pairs tested and their corresponding pvalues (using GLMM)
gaussian_glmm_qtl_results_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glmm_covariate_method_"$covariate_method"_permute_False_results.txt"
# All significant variant gene pairs (eFDR <= .01) and their corresponding pvalues (using GLMM)
gaussian_glmm_qtl_significant_results=$qtl_results_dir$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glmm_covariate_method_"$covariate_method"_efdr_.01_significant.txt"
# All significant genes (eFDR <= .01) and their top variant and their corresponding pvalues (using GLMM)
gaussian_glmm_qtl_significant_gene_results=$qtl_results_dir$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_glmm_covariate_method_"$covariate_method"_efdr_.01_significant_egenes.txt"


# Dynamic qtl results file from NB method
if [ $covariate_method == "none" ]; then
    nb_dynamic_qtl_file=$nb_dynamic_qtl_dir"te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_dosage_covariate_method_none_permutation_scheme_none_permute_False_merged_dynamic_qtl_results.txt"
else
    nb_dynamic_qtl_file=$nb_dynamic_qtl_dir"te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_dosage_covariate_method_cell_line_fixed_effect_"$covariate_method"Xtime_permutation_scheme_none_permute_False_merged_dynamic_qtl_results.txt"
fi


output_stem=$qtl_visualization_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_covariate_method_"$covariate_method
if false; then
sbatch temporal_qtl_enrichment_analysis.sh $gaussian_glm_qtl_results_file $gaussian_glm_qtl_significant_results $gaussian_glm_qtl_significant_gene_results $gaussian_glmm_qtl_results_file $gaussian_glmm_qtl_significant_results $gaussian_glmm_qtl_significant_gene_results $nb_dynamic_qtl_file $genotype_file $output_stem $t_0_eqtl_results_file
fi

if false; then
Rscript visualize_cell_line_overlap_between_methods.R $qtl_visualization_dir "200"
fi





##########################################
# Step 5: ChromHMM enrichment analysis
##########################################

num_permutations="1000"
model_version="glmm_time"
covariate_method="pc1_6"


if false; then
top_genes_arr=( "50" "75" "100" "125" "150" "175" "200" "300" "400" "500" "1000" "1500" "2000" "2500" "3000" "4000" "5000" "6000" )
for top_genes in "${top_genes_arr[@]}"
do
	parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
	real_dynamic_qtl_results_file=$qtl_results_dir$parameter_string"_permute_False_results.txt"
	echo $parameter_string
	sbatch chrom_hmm_enrichment_analysis_cross_top_genes.sh $parameter_string$top_genes"_" $num_permutations $chrom_hmm_input_dir $real_dynamic_qtl_results_file $time_step_independent_stem $chrom_hmm_enrichment_directory $top_genes
done
fi


base_parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version
if false; then
Rscript visualize_chrom_hmm_enrichment_cross_top_genes.R $base_parameter_string $num_permutations $chrom_hmm_enrichment_directory
fi






####################################################
# Step 6: Tissue specific chrom-hmm enrichments
####################################################
model_version="glm"
covariate_method="pc1_5"
num_permutations="1000"
num_genes="550"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method


variant_gene_pairs_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"


cluster_assignment_stem=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_sorted_variant_gene_pairs_"
##########################################################################
##########################################################################
cluster_assignment="dynamic_qtl_predicted_means_0.5"
##########################################################################
##########################################################################
if false; then
sh tissue_specific_chrom_hmm_enrichment_analysis.sh $parameter_string"_"$num_genes"_genes_" $num_permutations $chrom_hmm_input_dir $variant_gene_pairs_file $time_step_independent_stem $tissue_specific_chrom_hmm_enrichment_dir $cluster_assignment_stem $cluster_assignment $covariate_method $num_genes
fi


####################################################
# Step 7: Gene set enrichment analysis
####################################################
model_version="glm"
covariate_method="none"
num_genes="550"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method

variant_gene_pairs_file=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_genes_organized_input_data.txt"
cluster_assignment_stem=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_sorted_variant_gene_pairs_"
cluster_assignment="manual_annotations"
if false; then
sbatch gene_set_enrichment_analysis.sh $parameter_string $variant_gene_pairs_file $gencode_file $gene_set_enrichment_dir $gsea_data_dir $cluster_assignment_stem $cluster_assignment $time_step_independent_stem
fi



####################################################
# Step 8: Gene set enrichment analysis of cardiomyopathy genes
####################################################
model_version="glm"
covariate_method="none"
num_genes="550"

parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
variant_gene_pairs_file=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_genes_organized_input_data.txt"
all_hits_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"
if false; then
python cardiomyopathy_gene_set_enrichment_analysis.py $cardiomyopathy_gene_list $variant_gene_pairs_file $all_hits_file $gencode_file $cardiomyopathy_gene_set_enrichment_dir $parameter_string
fi



####################################################
# Step 9:GWAS Overlap
####################################################
model_version="glm_quadratic"
covariate_method="pc1_5"
num_genes="693"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
variant_gene_pairs_file=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_variant_gene_pairs_organized_input_data.txt"
all_hits_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"
if false; then
python compute_gwas_overlaps.py $gwas_catalog $gwas_overlap_dir $variant_gene_pairs_file $parameter_string $all_hits_file
fi


model_version="glm_quadratic"
covariate_method="pc1_5"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
num_genes="693"
parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method
variant_gene_pairs_file=$dynamic_qtl_visualization_dir$parameter_string"_top_"$num_genes"_variant_gene_pairs_organized_input_data.txt"
all_hits_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"
threshold="5e-8"
if false; then
python compute_gtex_gwas_overlaps.py $gtex_gwas_hits_dir $gwas_overlap_dir $variant_gene_pairs_file $parameter_string $all_hits_file $threshold
fi


####################################################
# Step 9: Time Step independent comparison file
####################################################
model_version="glm_quadratic"
covariate_method="pc1_5"
num_genes="693"
cluster_assignment="dynamic_qtl_quadratic_predicted_means_1.0"
# "dynamic_qtl_quadratic_predicted_means_1.0"
# "dynamic_qtl_predicted_means_1.0"


parameter_string="gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method


all_hits_file=$qtl_results_dir"gaussian_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_model_type_"$model_version"_covariate_method_"$covariate_method"_permute_False_results.txt"
if false; then
sh time_step_independent_comparison.sh $parameter_string"_top_"$num_genes"_" $time_step_independent_stem $time_step_comparison_dir $cluster_assignment $all_hits_file $num_genes
fi

if false; then
Rscript make_plots_for_manuscript.R $manuscript_figures_dir $dynamic_qtl_visualization_dir $time_step_comparison_dir $tissue_specific_chrom_hmm_enrichment_dir
fi




