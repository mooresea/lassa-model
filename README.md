# lassa-model
Estimate sub-national LASV infections and LF incidence rates in West Africa

This repository contains the code to recreate the modeling analysis in the article "Estimation of Lassa fever incidence rates in West Africa: development of a modeling framework to inform vaccine trial design". This article is currently available as a preprint at: https://doi.org/10.1101/2024.12.11.24318478 

The analysis involves running the following R scripts from the 'code' directory:<br/>
*Note: to successfully complete this analysis for the 2nd administrative level, the zip file containing population data, adm_2_pop_upd.z, will first need to be unzipped.<br />
<br/>
0. prior_iceberg.R<br />
    	Description: Estimate priors for proportion of infections that are asymptomatic/mild/unreported, a reported case, a reported death<br />
	Input: data/iceberg_data.csv<br />
	Output: results/prior_iceberg.RData<br />
2. step1_analyze_sero_lassa.R<br />
	Description: Estimate the FOI for each administrative unit with serology data using a catalytic model.<br />
	Input: data/case_reports_Lassa_IGG.csv<br />
	Output: results/foi_from_sero_adm_X_revrateXX.RData<br />
	Command variables:<br />
		1. admin: 1 or 2<br />
		2. reversion rate: 0,3,6 (%)<br />
3. step2_estimate_proportion_by_type_country_effects.R<br />
	Description: Estimate case and death reporting fractions for each country based on FOI estimates from step 1 and case/death data.<br />
	Input:  data/case_reports_Lassa_CASES.csv<br />
		results/foi_from_sero_admX_revrateXX.RData<br />
	Output: results/proportion_by_type_admX_revrateXX_country_upd.RData<br />
	Command variables:<br />
		1. admin: 1 or 2<br />
		2. reversion rate: 0,3,6 (%)<br />
4. step3_project_infections_by_country.R<br />
	Description: Estimate number of spillover infections for administrative units based on reporting fractions from step2.<br />
	Input:  data/case_reports_Lassa_CASES.csv<br />
			results/proportion_by_type_admX_revrateXX_country_upd.RData<br />
	Output: results/projected_infections_admX_revrateXX_country_upd.RData<br />
	Command variables:<br />
		1. admin: 1 or 2<br />
		2. reversion rate: 0,3,6 (%)<br />
5. step4_project_foi_by_country.R<br />
	Description: Project FOI for all adminstrative units based on estimated FOI rates.<br />
	Input:  results/foi_from_sero_adm_X_revrateXX.RData<br />
			results/projected_infections_admX_revrateXX_country_upd.RData<br />
	Output: results/projected_foi_admX_revrateXX_updated.RData<br />
	Command variables:<br />
		1. admin: 1 or 2<br />
		2. reversion rate: 0,3,6 (%)<br />
6. step5_simulate_spillovers_lassa_foi.R<br />
	Description: Project FOI for all adminstrative units based on estimated FOI rates.<br />
	Input:  results/projected_foi_admX_revrateXX_updated.RData<br />
	Output: results/simulated_spillovers_admX_revrateXX_',ifelse(raw.ind==1,'raw','modeled_outpredict'),'_serorev_risk.RData<br />
	Command variables:<br />
		1. admin: 1 or 2<br />
		2. reversion rate: 0,3,6 (%)<br />
		3. Indicator whether to use FOI estimates from step4 or from step4alt (ensemble model): 1 or 2<br />
7. step6_analyze_spillovers_lassa_foi.R<br />
	Description: Projection age-specific number of infections and LF cases for all adminstrative units for a set of reinfection risk scenarios. Output includes LF incidence rates reported in article Table 2 and supplemental tables.<br />
	Input: results/simulated_spillovers_admX_revrateXX_',ifelse(raw.ind==1,'raw','modeled_outpredict'),'_serorev_risk.RData<br />
	Output: results/LF_incidence_rates_admX_revrateXX_'ifelse(raw.ind==1,'_raw','_modeled_outpredict')'.csv<br />
	Command variables:<br />
		1. admin: 1 or 2<br />
		2. reversion rate: 0,3,6 (%)<br />
		3. Indicator whether to use FOI estimates from step4 or from step4alt (ensemble model): 1 or 2	<br />

<br />
To also complete the analysis that uses regression models to estimate FOI the following additional scripts should be run:<br />
A. step4alt_model_foi_estimation.R<br />
B. step4alt_model_covariate_significance.R<br />
C. step4alt_ensemble_weights_outpredict.R<br />
