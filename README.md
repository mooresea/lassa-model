# lassa-model
Estimate sub-national LASV infections and LF incidence rates in West Africa

This repository contains the code to recreate the modeling analysis in the article "Estimation of Lassa fever incidence rates in West Africa: development of a modeling framework to inform vaccine trial design". This article is currently available as a preprint at: https://doi.org/10.1101/2024.12.11.24318478 

The analysis involves running the following R scripts from the 'code' directory:

0. prior_iceberg.R
1. step1_analyze_sero_lassa.R
2. step2_estimate_proportion_by_type_country_effects.R
3. step3_project_infections_by_country.R
4. step4_project_foi_by_country.R
5. step5_simulate_spillovers_lassa_foi.R

*Note: to successfully complete this analysis for the 2nd administrative level, the zip file containing population data, adm_2_pop_upd.z, will first need to be unzipped.

To also complete the analysis that uses regression models to estimate FOI the following additional scripts should be run:
A. step4alt_model_foi_estimation.R
B. step4alt_model_covariate_significance.R
C. step4alt_ensemble_weights_outpredict.R
