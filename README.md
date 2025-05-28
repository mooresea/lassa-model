# lassa-model
Estimate sub-national LASV infections and LF incidence rates in West Africa

This repository contains the code to recreate the modeling analysis in the article "Estimation of Lassa fever incidence rates in West Africa: development of a modeling framework to inform vaccine trial design". This article is currently available as a preprint at: https://doi.org/10.1101/2024.12.11.24318478 

The analysis involves running the following R scripts from the 'code' directory:

0. prior_iceberg.R
	Description: Estimate priors for proportion of infections that are asymptomatic/mild/unreported, a reported case, a reported death
	Input: data/iceberg_data.csv
	Output: results/prior_iceberg.RData
1. step1_analyze_sero_lassa.R
	Description: Estimate the FOI for each administrative unit with serology data using a catalytic model.
	Input: data/case_reports_Lassa_IGG.csv
	Output: results/foi_from_sero_adm_X_revrateXX.RData
	Command variables:
		1. admin: 1 or 2
		2. reversion rate: 0,3,6 (%)
2. step2_estimate_proportion_by_type_country_effects.R
	Description: Estimate case and death reporting fractions for each country based on FOI estimates from step 1 and case/death data.
	Input:  data/case_reports_Lassa_CASES.csv
			results/foi_from_sero_admX_revrateXX.RData
	Output: results/proportion_by_type_admX_revrateXX_country_upd.RData
	Command variables:
		1. admin: 1 or 2
		2. reversion rate: 0,3,6 (%)
3. step3_project_infections_by_country.R
	Description: Estimate number of spillover infections for administrative units based on reporting fractions from step2.
	Input:  data/case_reports_Lassa_CASES.csv
			results/proportion_by_type_admX_revrateXX_country_upd.RData
	Output: results/projected_infections_admX_revrateXX_country_upd.RData
	Command variables:
		1. admin: 1 or 2
		2. reversion rate: 0,3,6 (%)
4. step4_project_foi_by_country.R
	Description: Project FOI for all adminstrative units based on estimated FOI rates.
	Input:  results/foi_from_sero_adm_X_revrateXX.RData
			results/projected_infections_admX_revrateXX_country_upd.RData
	Output: results/projected_foi_admX_revrateXX_updated.RData
	Command variables:
		1. admin: 1 or 2
		2. reversion rate: 0,3,6 (%)
5. step5_simulate_spillovers_lassa_foi.R
	Description: Project FOI for all adminstrative units based on estimated FOI rates.
	Input:  results/projected_foi_admX_revrateXX_updated.RData
	Output: results/simulated_spillovers_admX_revrateXX_',ifelse(raw.ind==1,'raw','modeled_outpredict'),'_serorev_risk.RData
	Command variables:
		1. admin: 1 or 2
		2. reversion rate: 0,3,6 (%)
		3. Indicator whether to use FOI estimates from step4 or from step4alt (ensemble model): 1 or 2
6. step6_analyze_spillovers_lassa_foi.R
	Description: Projection age-specific number of infections and LF cases for all adminstrative units for a set of reinfection risk scenarios. Output includes LF incidence rates reported in article Table 2 and supplemental tables.
	Input: results/simulated_spillovers_admX_revrateXX_',ifelse(raw.ind==1,'raw','modeled_outpredict'),'_serorev_risk.RData
	Output: results/LF_incidence_rates_admX_revrateXX_'ifelse(raw.ind==1,'_raw','_modeled_outpredict')'.csv
	Command variables:
		1. admin: 1 or 2
		2. reversion rate: 0,3,6 (%)
		3. Indicator whether to use FOI estimates from step4 or from step4alt (ensemble model): 1 or 2	
*Note: to successfully complete this analysis for the 2nd administrative level, the zip file containing population data, adm_2_pop_upd.z, will first need to be unzipped.

To also complete the analysis that uses regression models to estimate FOI the following additional scripts should be run:
A. step4alt_model_foi_estimation.R
B. step4alt_model_covariate_significance.R
C. step4alt_ensemble_weights_outpredict.R
