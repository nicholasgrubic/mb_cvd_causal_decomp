#################################################################################################################################################

# Project: Contribution of Adolescent Movement Behaviours to Socioeconomic Inequalities in Adults Cardiovascular Disease Risk
# Purpose: Function to perform causal decomposition for all interventional scenarios
# Author: Nicholas Grubic
# Creation Date: December 17, 2025
# Last Updated: April 27, 2026

#################################################################################################################################################

# Call packages
library(margins)
library(survey)
library(dplyr)
library(tidyr)
library(purrr)
library(GJRM)
library(mvtnorm)
library(msm)
library(svrep)

########################################################################################################################################
# IMPORT IMPUTED DATA
########################################################################################################################################

# load MI dataset
load("")

########################################################################################################################################
# INTERVENTIONAL SCENARIOS FUNCTION
# - Generates absolute (RD) and relative (RR) point estimates for natural course and counterfactual/intervention scenarios
# - Generates absolute (RD) and relative (RR) point estimates for observed data (for model misspecification check)
# - Also produces p0, p1, and population risk estimates for all scenarios and observed data
# - Function is applied to 1 imputed dataset (outputs list of matrices: point estimates for each BS replicate x MC draw)
# - A single natural course contrast is computed, along with joint and path-specific contrasts for both the elimination and equalization
#   scenarios
# - Function arguments:
#    1. data_imp = imputed dataset
#    2. n_mc = number of Monte-Carlo repetitions
#    3. bs = bootstrap toggle
#       * 'yes' if want point estimates for bootstrapped replicates
#       * 'no' if want point estimates for imputed dataset only (need to also specify n_boot = 1)
#    4. n_boot = number of bootstrap replicates
#    5. sep = SEP indicator (binary)
#    6. outcome = outcome (binary)
########################################################################################################################################

# Start of function
int_scenarios <- function(data_imp, n_mc, bs, n_boot, sep, outcome) {
  
  ######################################################################################################################################
  # PRE-ALLOCATION MATRICES
  ######################################################################################################################################
  
  # Create matrices to store RD, RR, p0, p1, and population risk estimates from nc and cf scenarios (joint and path-specific) for all Monte Carlo repetitions within each bootstrap replicate
  # Create vectors to store RD, RR, p0, p1, and population risk estimates from observed data for all bootstrap replicates
  
  # NATURAL COURSE
  RD_nc_coef_sim <- matrix(NA, nrow = n_boot, ncol = n_mc)
  RR_nc_coef_sim <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p0_nc_sim <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p1_nc_sim <- matrix(NA, nrow = n_boot, ncol = n_mc)
  pop_risk_nc_sim <- matrix(NA, nrow = n_boot, ncol = n_mc)
  
  # OBSERVED
  RD_ob_coef_sim <- rep(NA, n_boot)
  RR_ob_coef_sim <- rep(NA, n_boot)
  p0_ob_sim <- rep(NA, n_boot)
  p1_ob_sim <- rep(NA, n_boot)
  pop_risk_ob_sim <- rep(NA, n_boot)
  
  # ELIMINATION
  RD_cf_joint_coef_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  RR_cf_joint_coef_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p0_cf_joint_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p1_cf_joint_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  pop_risk_cf_joint_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  
  RD_cf_mvpa_coef_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  RR_cf_mvpa_coef_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p0_cf_mvpa_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p1_cf_mvpa_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  pop_risk_cf_mvpa_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  
  RD_cf_sleep_coef_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  RR_cf_sleep_coef_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p0_cf_sleep_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p1_cf_sleep_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  pop_risk_cf_sleep_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  
  RD_cf_st_coef_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  RR_cf_st_coef_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p0_cf_st_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p1_cf_st_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  pop_risk_cf_st_sim_el <- matrix(NA, nrow = n_boot, ncol = n_mc)
  
  # EQUALIZATION
  RD_cf_joint_coef_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  RR_cf_joint_coef_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p0_cf_joint_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p1_cf_joint_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  pop_risk_cf_joint_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  
  RD_cf_mvpa_coef_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  RR_cf_mvpa_coef_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p0_cf_mvpa_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p1_cf_mvpa_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  pop_risk_cf_mvpa_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  
  RD_cf_sleep_coef_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  RR_cf_sleep_coef_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p0_cf_sleep_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p1_cf_sleep_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  pop_risk_cf_sleep_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  
  RD_cf_st_coef_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  RR_cf_st_coef_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p0_cf_st_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  p1_cf_st_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  pop_risk_cf_st_sim_eq <- matrix(NA, nrow = n_boot, ncol = n_mc)
  
  ######################################################################################################################################
  # MI BOOT SET-UP - Schomaker & Heumann, 2019 (https://pmc.ncbi.nlm.nih.gov/articles/PMC5986623/)
  ######################################################################################################################################
  
  # STEPS:
  # - MI is performed (DONE)
  # - For each of the M imputed datasets, B bootstrap samples are drawn which yields B x M datasets
  # - In each of these data sets the measures of interest are estimated
  # - The pooled sample of ordered B x M estimates is used to construct the 95% confidence interval
  
  # Remove REGION + GSW12345 from imputed data (will be replaced from full GSW12345 respondents)
  link_i <- data_imp %>% dplyr::select(-REGION, -GSW12345)
  
  # Link imputed data back to full W12345 survey sample
  data_imp_link <- w12345_sample %>%
    dplyr::select(AID, REGION, PSUSCID, GSW12345, age_wave1) %>%
    left_join(link_i, by = "AID")
  
  # Build svydesign object to get bootstrap replicates
  bs_svydesign <- svydesign(
    ids = ~PSUSCID,
    strata = ~REGION,
    weights = ~GSW12345,
    nest = TRUE,
    survey.lonely.psu = "adjust",
    data = data_imp_link
  )
  
  # Generate bootstrap replicates on full W12345 survey sample
  # - Ensures replicate weights reflect PSU level variability, domain membership variability, and the original sampling design
  bootstrap_add.health <- as_bootstrap_design(bs_svydesign,
                                              type = "Rao-Wu-Yue-Beaumont",
                                              replicates = n_boot,
                                              compress = T,
                                              samp_method_by_stage = 'PPSWR')
  
  ######################################################################################################################################
  # START OF BOOTSTRAP
  # - After the bootstrap replicates are generated on the full W12345 survey sample, the weight for the kth replicate replaces GSW12345
  #   and the svydesign object for that replicate is generated (with appropriate subsetting to reflect the analytic sample)
  ######################################################################################################################################
  
  # Run bootstrap loop
  for (k in 1:n_boot) {
    
    if(bs == 'yes'){
      # Extract replicate weights for kth bootstrap
      GSW12345_bs <- bootstrap_add.health$repweights[,k] 
      
      # Extract survey variables
      svy_data_bs <- bootstrap_add.health$variables
      
      # Replace existing GSW12345 with the bootstrap replicate weight
      svy_data_bs$GSW12345 <- GSW12345_bs
      
      # Build svydesign object
      svydesign_imp <- svydesign(
        ids = ~PSUSCID,
        strata = ~REGION,
        weights = ~GSW12345,
        nest = TRUE,
        survey.lonely.psu = "adjust",
        data = svy_data_bs
      )
    }
    
    # If bs = 'no', the function outputs the point estimate (uses the original GSW12345 weight; also need to specify n_boot = 1)
    if(bs == 'no'){
      svydesign_imp <- bs_svydesign
    }
    
    # Apply analytic exclusions
    svydesign_imp <- subset(
      svydesign_imp,
      age_wave1 >= 12 &
        age_wave1 <= 17 &
        !is.na(age_wave1) &
        cvd_baseline == 'No' &
        # Zero-weight replicate weight clean-up: restrict to obs with GSW12345 > 0 (only relevant for svydesign with bootstrap replicate weights)
        # - Needed as svyolr breaks down when some obs have weight = 0 for some reason;
        # - Other svyglm models handle zero weights, so there is no difference in estimates/coefficients for these models whether or not this exclusion is added to the subset statement
        GSW12345 > 0
    )
    
    # Recode group (a) - need additional a_1, a_2, and a_3 for use in multivariate models 
    if(sep == 'par_edu'){
      decomp_svydesign <- update(
        svydesign_imp,
        a = ifelse(par_edu == 'College/University Degree or Higher', 0, 1),
        a_1 = ifelse(par_edu == 'College/University Degree or Higher', 0, 1),
        a_2 = ifelse(par_edu == 'College/University Degree or Higher', 0, 1),
        a_3 = ifelse(par_edu == 'College/University Degree or Higher', 0, 1),
        other_sep = fin_hardship) # For use in parametric models (additionally adjust for other SEP group)
    }
    
    if(sep == 'fin_hardship'){
      decomp_svydesign <- update(
        svydesign_imp,
        a = ifelse(fin_hardship == 'No', 0, 1),
        a_1 = ifelse(fin_hardship == 'No', 0, 1),
        a_2 = ifelse(fin_hardship == 'No', 0, 1),
        a_3 = ifelse(fin_hardship == 'No', 0, 1),
        other_sep = par_edu) # For use in parametric models (additionally adjust for other SEP group)
    }
    
    # Recode outcome (y)
    if(outcome == 'primary_outcome'){
      decomp_svydesign <- update(
        decomp_svydesign,
        y = ifelse(primary_outcome == 'Low-Intermediate Risk (<20%)', 0, 1))
    }
    
    if(outcome == 'hbp_w5'){
      decomp_svydesign <- update(
        decomp_svydesign,
        y = ifelse(hbp_w5 == 'No', 0, 1))
    }
    
    # Recode time-varying exposures at W1 and W3 (e1_x and e3_x) and time-varying covariates at W3 (l3_x)
    decomp_svydesign <- update(
      decomp_svydesign,
      e1_1 = ifelse(mvpa_w1 == 'Healthy', 0, 1),
      e1_2 = ifelse(sleep_w1 == 'Healthy', 0, 1),
      e1_3 = ifelse(screen_time_w1 == 'Healthy', 0, 1),
      e3_1 = ifelse(mvpa_w3 == 'Healthy', 0, 1),
      e3_2 = ifelse(sleep_w3 == 'Healthy', 0, 1),
      e3_3 = ifelse(screen_time_w3 == 'Healthy', 0, 1),
      l3_1 = ifelse(chronic_disease_w3 == 'No', 0, 1),
      l3_2 = ifelse(alcohol_misuse_w3 == 'No', 0, 1),
      l3_3 = ifelse(current_smoking_w3 == 'No', 0, 1),
      l3_4 = ifelse(mental_health_w3 == 'No', 0, 1),
      l3_5 = ifelse(obesity_w2 == 'No', 0, 1),
      l3_6 = factor(as.numeric(adj_hh_income_quart_w3), levels = 1:4))
    
    # Extract data from svydesign for predictions and simulations
    svy_data <- decomp_svydesign$variables
    
    # Create list of all stochastic variables and set to NA
    vars_to_na <- c("e1_1","e1_2","e1_3","e3_1","e3_2","e3_3","l3_1","l3_2","l3_3","l3_4","l3_5","l3_6")
    
    ######################################################################################################################################
    # G-FORMULA MODEL FITTING
    ######################################################################################################################################
    
    # Fit models for stochastic variables
    # - stochastic = outcome (y), time-varying exposures at W1 and W3 (e1_x and e3_x), and time-varying covariates at W3 (l3_x)
    # - deterministic/fixed (no arrows in DAG going into them) = a (group), c_x (baseline confounders), and l1_x (time-varying confounders at W1)
    # - wrap stochastic variables in factor() since when we simulate them, they generate a new variable that is considered numeric (simulated as 0 or 1); no difference for binary variables if no interaction in the model, but important if we include interactions
    # - all deterministic variables are already in proper numeric (continuous) or factor (categorical) form, so no need to wrap
    
    # Fit outcome model
    m_Y <- svyglm(y ~ 
                    factor(a) + #SEP group
                    factor(e1_1) + factor(e1_2) + factor(e1_3) + # exposures at W1
                    factor(e3_1) + factor(e3_2) + factor(e3_3) + # exposures at W3
                    age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                    chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1 + # all time-varying confounders at W1
                    factor(l3_1) + factor(l3_2) + factor(l3_3) + factor(l3_4) + factor(l3_5) + factor(l3_6), # all time-varying confounders at W3
                  family = quasibinomial(link = "logit"), 
                  design = decomp_svydesign)
    
    # Fit time-varying exposure models (W1 and W3) using multivariate models
    # W1 - trivariate probit model
    m_E1 <- gjrm(list(e1_1 ~ # mvpa at W1
                        factor(a_1) + # SEP group
                        age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                        chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1, # all time-varying confounders at W1
                      e1_2 ~ # sleep at W1
                        factor(a_2) + # SEP group
                        age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                        chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1, # all time-varying confounders at W1
                      e1_3 ~ # screen time at W1 
                        factor(a_3) + # SEP group
                        age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                        chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1), # all time-varying confounders at W1
                 data = svy_data,
                 weights = GSW12345,
                 margins = c("probit", "probit", "probit"),
                 model = "T")
    
    # Warning about "non-integer successes in a binomial glm" occurs because gjrm/binomial() interprets weights as trial counts
    # Our weights are survey/probability weights (not actual counts), so the warning can be safely ignored
    
    # Adjust covariance matrix from gjrm models fitted to complex survey data
    # - This function has been extracted from the survey package and adapted to the class of this package’s models
    # - It computes the sandwich variance estimator for a copula model fitted to the data from a complex sample survey (Lumley, 2004)
    
    # x = a fitted gjrm object as produced by the respective fitting function.
    # design = a svydesign object as produced by svydesign() from the survey package.
    m_E1 <- adjCovSD(x = m_E1, design = decomp_svydesign)
    
    # Extract estimated variance parameters and correlation coefficients of the error terms from the joint mediator model fit
    # Then, construct an estimated covariance matrix for the error terms using these components
    # For probit models, the latent variances are fixed at 1 (sigma1 = sigma2 = sigma3 = 1) and the correlations between error terms are given by the estimated theta values (see Smith et al., 2024 code)
    sigma1 <- 1
    sigma2 <- 1
    sigma3 <- 1
    
    rho12_m_E1 <- m_E1$theta12
    rho13_m_E1 <- m_E1$theta13
    rho23_m_E1 <- m_E1$theta23
    
    cov_mat_E1 <- matrix(c(
      sigma1,      rho12_m_E1,  rho13_m_E1,
      rho12_m_E1,  sigma2,      rho23_m_E1,
      rho13_m_E1,  rho23_m_E1,  sigma3
    ), nrow = 3, byrow = TRUE)
    
    # W3 - trivariate probit model
    m_E3 <- gjrm(list(e3_1 ~ # mvpa at W3
                        factor(a_1) + # SEP group
                        factor(e1_1) + factor(e1_2) + factor(e1_3) + # exposures at W1
                        age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                        chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1 + # all time-varying confounders at W1
                        factor(l3_1) + factor(l3_2) + factor(l3_3) + factor(l3_4) + factor(l3_5) + factor(l3_6), # all time-varying confounders at W3
                      e3_2 ~ # sleep at W3
                        factor(a_2) + # SEP group
                        factor(e1_1) + factor(e1_2) + factor(e1_3) + # exposures at W1
                        age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                        chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1 + # all time-varying confounders at W1
                        factor(l3_1) + factor(l3_2) + factor(l3_3) + factor(l3_4) + factor(l3_5) + factor(l3_6), # all time-varying confounders at W3
                      e3_3 ~  # screen time at W3
                        factor(a_3) + # SEP group
                        factor(e1_1) + factor(e1_2) + factor(e1_3) + # exposures at W1
                        age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                        chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1 + # all time-varying confounders at W1
                        factor(l3_1) + factor(l3_2) + factor(l3_3) + factor(l3_4) + factor(l3_5) + factor(l3_6)), # all time-varying confounders at W3
                 data = svy_data,
                 weights = GSW12345,
                 margins = c("probit", "probit", "probit"),
                 model = "T")
    
    m_E3 <- adjCovSD(x = m_E3, design = decomp_svydesign)
    
    sigma1 <- 1
    sigma2 <- 1
    sigma3 <- 1
    
    rho12_m_E3 <- m_E3$theta12
    rho13_m_E3 <- m_E3$theta13
    rho23_m_E3 <- m_E3$theta23
    
    cov_mat_E3 <- matrix(c(
      sigma1,      rho12_m_E3,  rho13_m_E3,
      rho12_m_E3,  sigma2,      rho23_m_E3,
      rho13_m_E3,  rho23_m_E3,  sigma3
    ), nrow = 3, byrow = TRUE)
    
    # Fit models for time-varying confounders at W3
    m_L3_1 <- svyglm(l3_1 ~ 
                       factor(a) + #SEP group
                       factor(e1_1) + factor(e1_2) + factor(e1_3) + # exposures at W1
                       age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                       chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1, # all prior time-varying confounders at W1
                     family = quasibinomial(link = "logit"),
                     design = decomp_svydesign)
    m_L3_2 <- svyglm(l3_2 ~ 
                       factor(a) + #SEP group
                       factor(e1_1) + factor(e1_2) + factor(e1_3) + # exposures at W1
                       age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                       chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1, # all prior time-varying confounders at W1
                     family = quasibinomial(link = "logit"),
                     design = decomp_svydesign)
    m_L3_3 <- svyglm(l3_3 ~ 
                       factor(a) + #SEP group
                       factor(e1_1) + factor(e1_2) + factor(e1_3) + # exposures at W1
                       age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                       chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1, # all prior time-varying confounders at W1
                     family = quasibinomial(link = "logit"),
                     design = decomp_svydesign)
    m_L3_4 <- svyglm(l3_4 ~ 
                       factor(a) + #SEP group
                       factor(e1_1) + factor(e1_2) + factor(e1_3) + # exposures at W1
                       age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                       chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1, # all prior time-varying confounders at W1
                     family = quasibinomial(link = "logit"),
                     design = decomp_svydesign)
    m_L3_5 <- svyglm(l3_5 ~ 
                       factor(a) + #SEP group
                       factor(e1_1) + factor(e1_2) + factor(e1_3) + # exposures at W1
                       age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                       chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1, # all prior time-varying confounders at W1
                     family = quasibinomial(link = "logit"),
                     design = decomp_svydesign)
    m_L3_6 <- svyolr(l3_6 ~ 
                       factor(a) + #SEP group
                       factor(e1_1) + factor(e1_2) + factor(e1_3) + # exposures at W1
                       age_group_wave1 + sex + race + adverse_birth + disability + popularity + adv_social_env + nhood_depr + parent_hhd_smoker + parent_alcohol_misuse + school_connect_score + family_connect_score + other_sep + # all baseline confounders
                       chronic_disease_w1 + alcohol_misuse_w1 + current_smoking_w1 + mental_health_w1 + obesity_w1 + adj_hh_income_quart_w1, # all prior time-varying confounders at W1)
                     design = decomp_svydesign)
    
    ######################################################################################################################################
    # SIMULATE NATURAL COURSE AND COUNTERFACTUAL/INTERVENTIONAL SCENARIOS
    ######################################################################################################################################
    
    # Monte Carlo details
    # - Each iteration of the Monte Carlo loop simulates all stochastic variables at all time points (in a sequential manner) and the downstream estimates are computed using that iteration of variables
    # - Repeat x times and average estimate over draws to get final estimates
    # - This approach preserves the joint distribution of variables conditional on the past (if separate loops per variable were performed, this would break the dependency structure between new variables and past variables, which could bias the estimate)
    # - Run SEPARATE loops for the natural course and each counterfactual/interventional scenarios (as done by Smith et al. 2025)
    # - The nc stochastic draws should not be influenced by the cf stochastic draws or vice-versa (can sometimes happen if both scenarios are in the same Monte Carlo loop, although it is minimized when # of MC draws are sufficiently large)
    
    #### COMPUTE NATURAL COURSE OUTCOME (y_nc) ###
    
    for (i in 1:n_mc) {
      
      # Create dataset for natural course scenario
      d_nc <- svy_data %>% mutate(across(all_of(vars_to_na), ~ NA))
      
      # Predict probabilities of e1_x (exposure W1) and simulate natural course of e1_x
      # Generate error terms for each group's natural outcomes
      eps_nc_E1 <- rmvnorm(nrow(d_nc), mean = c(0,0,0), sigma = cov_mat_E1)
      d_nc$e1_1 <- ifelse(predict(m_E1, newdata = d_nc, eq = 1) + eps_nc_E1[,1] > 0, 1, 0)
      d_nc$e1_2 <- ifelse(predict(m_E1, newdata = d_nc, eq = 2) + eps_nc_E1[,2] > 0, 1, 0)
      d_nc$e1_3 <- ifelse(predict(m_E1, newdata = d_nc, eq = 3) + eps_nc_E1[,3] > 0, 1, 0)
      
      # Predict probability of l3_x (time-varying confounders W3) with simulated e1_x
      d_nc$l3_1 <- rbinom(nrow(d_nc), 1, predict(m_L3_1, newdata = d_nc, type = "response"))
      d_nc$l3_2 <- rbinom(nrow(d_nc), 1, predict(m_L3_2, newdata = d_nc, type = "response"))
      d_nc$l3_3 <- rbinom(nrow(d_nc), 1, predict(m_L3_3, newdata = d_nc, type = "response"))
      d_nc$l3_4 <- rbinom(nrow(d_nc), 1, predict(m_L3_4, newdata = d_nc, type = "response"))
      d_nc$l3_5 <- rbinom(nrow(d_nc), 1, predict(m_L3_5, newdata = d_nc, type = "response"))
      
      p_mat <- predict(m_L3_6, newdata = d_nc, type = "probs")
      d_nc$l3_6 <- apply(p_mat, 1, function(p) {
        sample(1:4, size = 1, prob = p)})
      d_nc$l3_6 <- as.factor(d_nc$l3_6)
      
      # Predict probability of e3_x (exposure W3) with simulated e1_x and l3_x and simulate natural course of e3_x
      # Generate error terms for each group's natural outcomes
      eps_nc_E3 <- rmvnorm(nrow(d_nc), mean = c(0,0,0), sigma = cov_mat_E3)
      d_nc$e3_1 <- ifelse(predict(m_E3, newdata = d_nc, eq = 1) + eps_nc_E3[,1] > 0, 1, 0)
      d_nc$e3_2 <- ifelse(predict(m_E3, newdata = d_nc, eq = 2) + eps_nc_E3[,2] > 0, 1, 0)
      d_nc$e3_3 <- ifelse(predict(m_E3, newdata = d_nc, eq = 3) + eps_nc_E3[,3] > 0, 1, 0)
      
      # Predict probability of y (outcome) with simulated e1_x, e3_x, and l3_x and simulate natural course of y
      d_nc$y_nc <- rbinom(nrow(d_nc), 1, predict(m_Y, newdata = d_nc, type = "response"))
      
      # Add simulated nc outcome to svydesign object
      decomp_svydesign_nc <- update(
        decomp_svydesign,
        y_nc = d_nc$y_nc)
      
      # Compute CDM (RD and RR) for natural course scenario = use y_nc
      final_model_nc <- svyglm(y_nc ~ age_group_wave1 + sex, family = quasibinomial(link = "logit"), design = decomp_svydesign_nc)
      risks_nc <- svypredmeans(adjustmodel = final_model_nc, groupfactor = ~ a)
      RD_nc_contrast <- svycontrast(risks_nc, quote(`1` - `0`))
      RR_nc_contrast <- svycontrast(risks_nc, quote(`1` / `0`))
      RD_nc_coef_sim[k, i] <- as.numeric(coef(RD_nc_contrast))
      RR_nc_coef_sim[k, i] <- as.numeric(coef(RR_nc_contrast))
      
      # Compute group-specific and population-wide risks
      p0_nc_sim[k, i] <- as.numeric(coef(risks_nc)[1])
      p1_nc_sim[k, i] <- as.numeric(coef(risks_nc)[2])
      
      pred_probs_nc <- predict(final_model_nc, type = "response")  
      pop_risk_nc <- svymean(pred_probs_nc, design = decomp_svydesign_nc)
      pop_risk_nc_sim[k, i] <- as.numeric(coef(pop_risk_nc)[1])
    }
    
    #### COMPUTE COUNTERFACTUAL/INTERVENTION OUTCOMES (y_cf) ###
    
    #### ELIMINATION ###
    
    # JOINT
    
    for (i in 1:n_mc) {  
      
      # Create dataset for counterfactual/intervention scenario
      d_cf_joint <- svy_data %>% mutate(across(all_of(vars_to_na), ~ NA))
      
      # Set e1_x and e3_x - determined from joint intervention (everyone is healthy at both time points)
      d_cf_joint$e1_1 <- 0
      d_cf_joint$e1_2 <- 0
      d_cf_joint$e1_3 <- 0
      d_cf_joint$e3_1 <- 0
      d_cf_joint$e3_2 <- 0
      d_cf_joint$e3_3 <- 0
      
      # Predict probability of l3_x (time-varying confounders w3) with assigned e1_x
      d_cf_joint$l3_1 <- rbinom(nrow(d_cf_joint), 1, predict(m_L3_1, newdata = d_cf_joint, type = "response"))
      d_cf_joint$l3_2 <- rbinom(nrow(d_cf_joint), 1, predict(m_L3_2, newdata = d_cf_joint, type = "response"))
      d_cf_joint$l3_3 <- rbinom(nrow(d_cf_joint), 1, predict(m_L3_3, newdata = d_cf_joint, type = "response"))
      d_cf_joint$l3_4 <- rbinom(nrow(d_cf_joint), 1, predict(m_L3_4, newdata = d_cf_joint, type = "response"))
      d_cf_joint$l3_5 <- rbinom(nrow(d_cf_joint), 1, predict(m_L3_5, newdata = d_cf_joint, type = "response"))
      
      p_mat <- predict(m_L3_6, newdata = d_cf_joint, type = "probs")
      d_cf_joint$l3_6 <- apply(p_mat, 1, function(p) {
        sample(1:4, size = 1, prob = p)})
      d_cf_joint$l3_6 <- as.factor(d_cf_joint$l3_6)
      
      # Predict probability of y (outcome) with assigned e1_x and e3_x, and simulated l3_x and simulate natural course of y
      d_cf_joint$y_cf_joint <- rbinom(nrow(d_cf_joint), 1, predict(m_Y, newdata = d_cf_joint, type = "response"))
      
      # Add simulated cf scenario outcome to svydesign object
      decomp_svydesign_cf_joint <- update(
        decomp_svydesign,
        y_cf_joint = d_cf_joint$y_cf_joint
      )
      
      # Compute CDM (RD and RR) for cf scenario = use y_cf
      final_model_cf_joint <- svyglm(y_cf_joint ~ age_group_wave1 + sex, family = quasibinomial(link = "logit"), design = decomp_svydesign_cf_joint)
      risks_cf_joint <- svypredmeans(adjustmodel = final_model_cf_joint, groupfactor = ~ a)
      RD_cf_joint_contrast <- svycontrast(risks_cf_joint, quote(`1` - `0`))
      RR_cf_joint_contrast <- svycontrast(risks_cf_joint, quote(`1` / `0`))
      RD_cf_joint_coef_sim_el[k, i] <- as.numeric(coef(RD_cf_joint_contrast)) 
      RR_cf_joint_coef_sim_el[k, i] <- as.numeric(coef(RR_cf_joint_contrast))
      
      # Compute group-specific and population-wide risks
      p0_cf_joint_sim_el[k, i] <- as.numeric(coef(risks_cf_joint)[1])
      p1_cf_joint_sim_el[k, i] <- as.numeric(coef(risks_cf_joint)[2])
      
      pred_probs_cf_joint <- predict(final_model_cf_joint, type = "response")  
      pop_risk_cf_joint <- svymean(pred_probs_cf_joint, design = decomp_svydesign_cf_joint)
      pop_risk_cf_joint_sim_el[k, i] <- as.numeric(coef(pop_risk_cf_joint)[1])
    }
    
    # PATH-SPECIFIC (MVPA)
    
    for (i in 1:n_mc) {  
      
      # Create dataset for counterfactual/intervention scenario
      d_cf_mvpa <- svy_data %>% mutate(across(all_of(vars_to_na), ~ NA))
      
      # Set e1_x and e3_x - determined from path-specific intervention (everyone is healthy for specific MB at both time points)
      d_cf_mvpa$e1_1 <- 0
      d_cf_mvpa$e3_1 <- 0
      
      # Predict probabilities of other e1_x (exposure W1) and simulate their natural courses
      # Generate error terms for each group's natural outcomes
      eps_cf_mvpa_E1 <- rmvnorm(nrow(d_cf_mvpa), mean = c(0,0,0), sigma = cov_mat_E1)
      d_cf_mvpa$e1_2 <- ifelse(predict(m_E1, newdata = d_cf_mvpa, eq = 2) + eps_cf_mvpa_E1[,2] > 0, 1, 0)
      d_cf_mvpa$e1_3 <- ifelse(predict(m_E1, newdata = d_cf_mvpa, eq = 3) + eps_cf_mvpa_E1[,3] > 0, 1, 0)
      
      # Predict probability of l3_x (time-varying confounders w3) with assigned e1_x
      d_cf_mvpa$l3_1 <- rbinom(nrow(d_cf_mvpa), 1, predict(m_L3_1, newdata = d_cf_mvpa, type = "response"))
      d_cf_mvpa$l3_2 <- rbinom(nrow(d_cf_mvpa), 1, predict(m_L3_2, newdata = d_cf_mvpa, type = "response"))
      d_cf_mvpa$l3_3 <- rbinom(nrow(d_cf_mvpa), 1, predict(m_L3_3, newdata = d_cf_mvpa, type = "response"))
      d_cf_mvpa$l3_4 <- rbinom(nrow(d_cf_mvpa), 1, predict(m_L3_4, newdata = d_cf_mvpa, type = "response"))
      d_cf_mvpa$l3_5 <- rbinom(nrow(d_cf_mvpa), 1, predict(m_L3_5, newdata = d_cf_mvpa, type = "response"))
      
      p_mat <- predict(m_L3_6, newdata = d_cf_mvpa, type = "probs")
      d_cf_mvpa$l3_6 <- apply(p_mat, 1, function(p) {
        sample(1:4, size = 1, prob = p)})
      d_cf_mvpa$l3_6 <- as.factor(d_cf_mvpa$l3_6)
      
      # Predict probabilities of other e3_x (exposure w3) and simulate their natural courses (path-specific only)
      # Generate error terms for each group's natural outcomes
      eps_cf_mvpa_E3 <- rmvnorm(nrow(d_cf_mvpa), mean = c(0,0,0), sigma = cov_mat_E3)
      d_cf_mvpa$e3_2 <- ifelse(predict(m_E3, newdata = d_cf_mvpa, eq = 2) + eps_cf_mvpa_E3[,2] > 0, 1, 0)
      d_cf_mvpa$e3_3 <- ifelse(predict(m_E3, newdata = d_cf_mvpa, eq = 3) + eps_cf_mvpa_E3[,3] > 0, 1, 0)
      
      # Predict probability of y (outcome) with assigned or simulated e1_x and e3_x, and simulated l3_x and simulate natural course of y
      d_cf_mvpa$y_cf_mvpa <- rbinom(nrow(d_cf_mvpa), 1, predict(m_Y, newdata = d_cf_mvpa, type = "response"))
      
      # Add simulated cf scenario outcomes to svydesign object
      decomp_svydesign_cf_mvpa <- update(
        decomp_svydesign,
        y_cf_mvpa = d_cf_mvpa$y_cf_mvpa
      )
      
      # Compute CDM (RD and RR) for cf scenario = use y_cf
      final_model_cf_mvpa <- svyglm(y_cf_mvpa ~ age_group_wave1 + sex, family = quasibinomial(link = "logit"), design = decomp_svydesign_cf_mvpa)
      risks_cf_mvpa <- svypredmeans(adjustmodel = final_model_cf_mvpa, groupfactor = ~ a)
      RD_cf_mvpa_contrast <- svycontrast(risks_cf_mvpa, quote(`1` - `0`))
      RR_cf_mvpa_contrast <- svycontrast(risks_cf_mvpa, quote(`1` / `0`))
      RD_cf_mvpa_coef_sim_el[k, i] <- as.numeric(coef(RD_cf_mvpa_contrast)) 
      RR_cf_mvpa_coef_sim_el[k, i] <- as.numeric(coef(RR_cf_mvpa_contrast)) 
      
      # Compute group-specific and population-wide risks
      p0_cf_mvpa_sim_el[k, i] <- as.numeric(coef(risks_cf_mvpa)[1])
      p1_cf_mvpa_sim_el[k, i] <- as.numeric(coef(risks_cf_mvpa)[2])
      
      pred_probs_cf_mvpa <- predict(final_model_cf_mvpa, type = "response")  
      pop_risk_cf_mvpa <- svymean(pred_probs_cf_mvpa, design = decomp_svydesign_cf_mvpa)
      pop_risk_cf_mvpa_sim_el[k, i] <- as.numeric(coef(pop_risk_cf_mvpa)[1])
    }
    
    # PATH-SPECIFIC (SLEEP)
    
    for (i in 1:n_mc) {  
      
      # Create dataset for counterfactual/intervention scenario
      d_cf_sleep <- svy_data %>% mutate(across(all_of(vars_to_na), ~ NA))
      
      # Set e1_x and e3_x - determined from path-specific intervention (everyone is healthy for specific MB at both time points)
      d_cf_sleep$e1_2 <- 0
      d_cf_sleep$e3_2 <- 0
      
      # Predict probabilities of other e1_x (exposure W1) and simulate their natural courses
      # Generate error terms for each group's natural outcomes
      eps_cf_sleep_E1 <- rmvnorm(nrow(d_cf_sleep), mean = c(0,0,0), sigma = cov_mat_E1)
      d_cf_sleep$e1_1 <- ifelse(predict(m_E1, newdata = d_cf_sleep, eq = 1) + eps_cf_sleep_E1[,1] > 0, 1, 0)
      d_cf_sleep$e1_3 <- ifelse(predict(m_E1, newdata = d_cf_sleep, eq = 3) + eps_cf_sleep_E1[,3] > 0, 1, 0)
      
      # Predict probability of l3_x (time-varying confounders w3) with assigned e1_x
      d_cf_sleep$l3_1 <- rbinom(nrow(d_cf_sleep), 1, predict(m_L3_1, newdata = d_cf_sleep, type = "response"))
      d_cf_sleep$l3_2 <- rbinom(nrow(d_cf_sleep), 1, predict(m_L3_2, newdata = d_cf_sleep, type = "response"))
      d_cf_sleep$l3_3 <- rbinom(nrow(d_cf_sleep), 1, predict(m_L3_3, newdata = d_cf_sleep, type = "response"))
      d_cf_sleep$l3_4 <- rbinom(nrow(d_cf_sleep), 1, predict(m_L3_4, newdata = d_cf_sleep, type = "response"))
      d_cf_sleep$l3_5 <- rbinom(nrow(d_cf_sleep), 1, predict(m_L3_5, newdata = d_cf_sleep, type = "response"))
      
      p_mat <- predict(m_L3_6, newdata = d_cf_sleep, type = "probs")
      d_cf_sleep$l3_6 <- apply(p_mat, 1, function(p) {
        sample(1:4, size = 1, prob = p)})
      d_cf_sleep$l3_6 <- as.factor(d_cf_sleep$l3_6)
      
      # Predict probabilities of other e3_x (exposure w3) and simulate their natural courses (path-specific only)
      # Generate error terms for each group's natural outcomes
      eps_cf_sleep_E3 <- rmvnorm(nrow(d_cf_sleep), mean = c(0,0,0), sigma = cov_mat_E3)
      d_cf_sleep$e3_1 <- ifelse(predict(m_E3, newdata = d_cf_sleep, eq = 1) + eps_cf_sleep_E3[,1] > 0, 1, 0)
      d_cf_sleep$e3_3 <- ifelse(predict(m_E3, newdata = d_cf_sleep, eq = 3) + eps_cf_sleep_E3[,3] > 0, 1, 0)
      
      # Predict probability of y (outcome) with assigned or simulated e1_x and e3_x, and simulated l3_x and simulate natural course of y
      d_cf_sleep$y_cf_sleep <- rbinom(nrow(d_cf_sleep), 1, predict(m_Y, newdata = d_cf_sleep, type = "response"))
      
      # Add simulated cf scenario outcomes to svydesign object
      decomp_svydesign_cf_sleep <- update(
        decomp_svydesign,
        y_cf_sleep = d_cf_sleep$y_cf_sleep
      )
      
      # Compute CDM (RD and RR) for cf scenario = use y_cf
      final_model_cf_sleep <- svyglm(y_cf_sleep ~ age_group_wave1 + sex, family = quasibinomial(link = "logit"), design = decomp_svydesign_cf_sleep)
      risks_cf_sleep <- svypredmeans(adjustmodel = final_model_cf_sleep, groupfactor = ~ a)
      RD_cf_sleep_contrast <- svycontrast(risks_cf_sleep, quote(`1` - `0`))
      RR_cf_sleep_contrast <- svycontrast(risks_cf_sleep, quote(`1` / `0`))
      RD_cf_sleep_coef_sim_el[k, i] <- as.numeric(coef(RD_cf_sleep_contrast)) 
      RR_cf_sleep_coef_sim_el[k, i] <- as.numeric(coef(RR_cf_sleep_contrast)) 
      
      # Compute group-specific and population-wide risks
      p0_cf_sleep_sim_el[k, i] <- as.numeric(coef(risks_cf_sleep)[1])
      p1_cf_sleep_sim_el[k, i] <- as.numeric(coef(risks_cf_sleep)[2])
      
      pred_probs_cf_sleep <- predict(final_model_cf_sleep, type = "response")  
      pop_risk_cf_sleep <- svymean(pred_probs_cf_sleep, design = decomp_svydesign_cf_sleep)
      pop_risk_cf_sleep_sim_el[k, i] <- as.numeric(coef(pop_risk_cf_sleep)[1])
    }
    
    # PATH-SPECIFIC (SCREEN TIME)
    
    for (i in 1:n_mc) {
      
      # Create dataset for counterfactual/intervention scenario
      d_cf_st <- svy_data %>% mutate(across(all_of(vars_to_na), ~ NA))
      
      # Set e1_x and e3_x - determined from path-specific intervention (everyone is healthy for specific MB at both time points)
      d_cf_st$e1_3 <- 0
      d_cf_st$e3_3 <- 0
      
      # Predict probabilities of other e1_x (exposure W1) and simulate their natural courses
      # Generate error terms for each group's natural outcomes
      eps_cf_st_E1 <- rmvnorm(nrow(d_cf_st), mean = c(0,0,0), sigma = cov_mat_E1)
      d_cf_st$e1_1 <- ifelse(predict(m_E1, newdata = d_cf_st, eq = 1) + eps_cf_st_E1[,1] > 0, 1, 0)
      d_cf_st$e1_2 <- ifelse(predict(m_E1, newdata = d_cf_st, eq = 2) + eps_cf_st_E1[,2] > 0, 1, 0)
      
      # Predict probability of l3_x (time-varying confounders w3) with assigned e1_x
      d_cf_st$l3_1 <- rbinom(nrow(d_cf_st), 1, predict(m_L3_1, newdata = d_cf_st, type = "response"))
      d_cf_st$l3_2 <- rbinom(nrow(d_cf_st), 1, predict(m_L3_2, newdata = d_cf_st, type = "response"))
      d_cf_st$l3_3 <- rbinom(nrow(d_cf_st), 1, predict(m_L3_3, newdata = d_cf_st, type = "response"))
      d_cf_st$l3_4 <- rbinom(nrow(d_cf_st), 1, predict(m_L3_4, newdata = d_cf_st, type = "response"))
      d_cf_st$l3_5 <- rbinom(nrow(d_cf_st), 1, predict(m_L3_5, newdata = d_cf_st, type = "response"))
      
      p_mat <- predict(m_L3_6, newdata = d_cf_st, type = "probs")
      d_cf_st$l3_6 <- apply(p_mat, 1, function(p) {
        sample(1:4, size = 1, prob = p)})
      d_cf_st$l3_6 <- as.factor(d_cf_st$l3_6)
      
      # Predict probabilities of other e3_x (exposure w3) and simulate their natural courses (path-specific only)
      # Generate error terms for each group's natural outcomes
      eps_cf_st_E3 <- rmvnorm(nrow(d_cf_st), mean = c(0,0,0), sigma = cov_mat_E3)
      d_cf_st$e3_1 <- ifelse(predict(m_E3, newdata = d_cf_st, eq = 1) + eps_cf_st_E3[,1] > 0, 1, 0)
      d_cf_st$e3_2 <- ifelse(predict(m_E3, newdata = d_cf_st, eq = 2) + eps_cf_st_E3[,2] > 0, 1, 0)
      
      # Predict probability of y (outcome) with assigned or simulated e1_x and e3_x, and simulated l3_x and simulate natural course of y
      d_cf_st$y_cf_st <- rbinom(nrow(d_cf_st), 1, predict(m_Y, newdata = d_cf_st, type = "response"))
      
      # Add simulated cf scenario outcomes to svydesign object
      decomp_svydesign_cf_st <- update(
        decomp_svydesign,
        y_cf_st = d_cf_st$y_cf_st
      )
      
      # Compute CDM (RD and RR) for cf scenario = use y_cf
      final_model_cf_st <- svyglm(y_cf_st ~ age_group_wave1 + sex, family = quasibinomial(link = "logit"), design = decomp_svydesign_cf_st)
      risks_cf_st <- svypredmeans(adjustmodel = final_model_cf_st, groupfactor = ~ a)
      RD_cf_st_contrast <- svycontrast(risks_cf_st, quote(`1` - `0`))
      RR_cf_st_contrast <- svycontrast(risks_cf_st, quote(`1` / `0`))
      RD_cf_st_coef_sim_el[k, i] <- as.numeric(coef(RD_cf_st_contrast)) 
      RR_cf_st_coef_sim_el[k, i] <- as.numeric(coef(RR_cf_st_contrast)) 
      
      # Compute group-specific and population-wide risks
      p0_cf_st_sim_el[k, i] <- as.numeric(coef(risks_cf_st)[1])
      p1_cf_st_sim_el[k, i] <- as.numeric(coef(risks_cf_st)[2])
      
      pred_probs_cf_st <- predict(final_model_cf_st, type = "response")  
      pop_risk_cf_st <- svymean(pred_probs_cf_st, design = decomp_svydesign_cf_st)
      pop_risk_cf_st_sim_el[k, i] <- as.numeric(coef(pop_risk_cf_st)[1])
    }
    
    ### EQUALIZATION ###
    
    # JOINT
    
    for (i in 1:n_mc) {  
      
      # Create dataset for counterfactual/intervention scenario
      d_cf_joint <- svy_data %>% mutate(across(all_of(vars_to_na), ~ NA))
      
      # Set e1_x - determined from joint intervention (distributions of all behaviours in a=1 group are equalized to those in a=0 group)
      # Set entire population to unexposed (a1=0, a2=0, a3=0, needed for multivariate exposure model)
      d_cf_joint_e1 <- d_cf_joint
      d_cf_joint_e1[c("a_1", "a_2", "a_3")] <- 0
      
      # Predict probabilities of e1_x (exposure W1) and simulate their values
      # Generate error terms for each group's counterfactual outcomes
      eps_cf_joint_E1 <- rmvnorm(nrow(d_cf_joint_e1), mean = c(0,0,0), sigma = cov_mat_E1)
      d_cf_joint_e1$e1_1 <- ifelse(predict(m_E1, newdata = d_cf_joint_e1, eq = 1) + eps_cf_joint_E1[,1] > 0, 1, 0)
      d_cf_joint_e1$e1_2 <- ifelse(predict(m_E1, newdata = d_cf_joint_e1, eq = 2) + eps_cf_joint_E1[,2] > 0, 1, 0)
      d_cf_joint_e1$e1_3 <- ifelse(predict(m_E1, newdata = d_cf_joint_e1, eq = 3) + eps_cf_joint_E1[,3] > 0, 1, 0)
      
      # Add simulated e1_x to dataset for future predictions + simulations
      d_cf_joint[c("e1_1", "e1_2", "e1_3")] <- d_cf_joint_e1[c("e1_1", "e1_2", "e1_3")]
      
      # Predict probability of l3_x (time-varying confounders w3) with assigned e1_x
      d_cf_joint$l3_1 <- rbinom(nrow(d_cf_joint), 1, predict(m_L3_1, newdata = d_cf_joint, type = "response"))
      d_cf_joint$l3_2 <- rbinom(nrow(d_cf_joint), 1, predict(m_L3_2, newdata = d_cf_joint, type = "response"))
      d_cf_joint$l3_3 <- rbinom(nrow(d_cf_joint), 1, predict(m_L3_3, newdata = d_cf_joint, type = "response"))
      d_cf_joint$l3_4 <- rbinom(nrow(d_cf_joint), 1, predict(m_L3_4, newdata = d_cf_joint, type = "response"))
      d_cf_joint$l3_5 <- rbinom(nrow(d_cf_joint), 1, predict(m_L3_5, newdata = d_cf_joint, type = "response"))
      
      p_mat <- predict(m_L3_6, newdata = d_cf_joint, type = "probs")
      d_cf_joint$l3_6 <- apply(p_mat, 1, function(p) {
        sample(1:4, size = 1, prob = p)})
      d_cf_joint$l3_6 <- as.factor(d_cf_joint$l3_6)
      
      # Set e3_x - determined from joint intervention (distributions of all behaviours in a=1 group are equalized to those in a=0 group)
      # Set entire population to unexposed (a1=0, a2=0, a3=0, needed for multivariate exposure model)
      d_cf_joint_e3 <- d_cf_joint
      d_cf_joint_e3[c("a_1", "a_2", "a_3")] <- 0
      
      # Predict probabilities of e3_x (exposure W3) and simulate their values
      # Generate error terms for each group's counterfactual outcomes
      eps_cf_joint_E3 <- rmvnorm(nrow(d_cf_joint_e3), mean = c(0,0,0), sigma = cov_mat_E3)
      d_cf_joint_e3$e3_1 <- ifelse(predict(m_E3, newdata = d_cf_joint_e3, eq = 1) + eps_cf_joint_E3[,1] > 0, 1, 0)
      d_cf_joint_e3$e3_2 <- ifelse(predict(m_E3, newdata = d_cf_joint_e3, eq = 2) + eps_cf_joint_E3[,2] > 0, 1, 0)
      d_cf_joint_e3$e3_3 <- ifelse(predict(m_E3, newdata = d_cf_joint_e3, eq = 3) + eps_cf_joint_E3[,3] > 0, 1, 0)
      
      # Add simulated e3_x to dataset for future predictions + simulations
      d_cf_joint[c("e3_1", "e3_2", "e3_3")] <- d_cf_joint_e3[c("e3_1", "e3_2", "e3_3")]
      
      # Predict probability of y (outcome) with simulated e1_x and e3_x, and simulated l3_x and simulate y
      d_cf_joint$y_cf_joint <- rbinom(nrow(d_cf_joint), 1, predict(m_Y, newdata = d_cf_joint, type = "response"))
      
      # Add simulated cf scenario outcome to svydesign object
      decomp_svydesign_cf_joint <- update(
        decomp_svydesign,
        y_cf_joint = d_cf_joint$y_cf_joint
      )
      
      # Compute CDM (RD and RR) for cf scenario = use y_cf
      final_model_cf_joint <- svyglm(y_cf_joint ~ age_group_wave1 + sex, family = quasibinomial(link = "logit"), design = decomp_svydesign_cf_joint)
      risks_cf_joint <- svypredmeans(adjustmodel = final_model_cf_joint, groupfactor = ~ a)
      RD_cf_joint_contrast <- svycontrast(risks_cf_joint, quote(`1` - `0`))
      RR_cf_joint_contrast <- svycontrast(risks_cf_joint, quote(`1` / `0`))
      RD_cf_joint_coef_sim_eq[k, i] <- as.numeric(coef(RD_cf_joint_contrast))
      RR_cf_joint_coef_sim_eq[k, i] <- as.numeric(coef(RR_cf_joint_contrast)) 
      
      # Compute group-specific and population-wide risks
      p0_cf_joint_sim_eq[k, i] <- as.numeric(coef(risks_cf_joint)[1])
      p1_cf_joint_sim_eq[k, i] <- as.numeric(coef(risks_cf_joint)[2])
      
      pred_probs_cf_joint <- predict(final_model_cf_joint, type = "response")  
      pop_risk_cf_joint <- svymean(pred_probs_cf_joint, design = decomp_svydesign_cf_joint)
      pop_risk_cf_joint_sim_eq[k, i] <- as.numeric(coef(pop_risk_cf_joint)[1])
    }
    
    # PATH-SPECIFIC (MVPA)
    
    for (i in 1:n_mc) {  
      
      # Create dataset for counterfactual/intervention scenario
      d_cf_mvpa <- svy_data %>% mutate(across(all_of(vars_to_na), ~ NA))
      
      # Set e1_x - determined from path-specific intervention (distribution of specific MB in a=1 group is equalized to those in a=0 group)
      # Set entire population to unexposed for specific MB (a1=0 [MVPA], a2=1 [screen time], a3=1 [sleep], needed for multivariate exposure model)
      d_cf_mvpa_e1 <- d_cf_mvpa
      d_cf_mvpa_e1$a_1 <- 0
      
      # Predict probability of e1_x (exposure W1) and simulate value
      # Generate error terms for each group's natural/counterfactual outcomes
      eps_cf_mvpa_E1 <- rmvnorm(nrow(d_cf_mvpa_e1), mean = c(0,0,0), sigma = cov_mat_E1)
      d_cf_mvpa_e1$e1_1 <- ifelse(predict(m_E1, newdata = d_cf_mvpa_e1, eq = 1) + eps_cf_mvpa_E1[,1] > 0, 1, 0)
      d_cf_mvpa_e1$e1_2 <- ifelse(predict(m_E1, newdata = d_cf_mvpa_e1, eq = 2) + eps_cf_mvpa_E1[,2] > 0, 1, 0)
      d_cf_mvpa_e1$e1_3 <- ifelse(predict(m_E1, newdata = d_cf_mvpa_e1, eq = 3) + eps_cf_mvpa_E1[,3] > 0, 1, 0)
      
      # Add simulated e1_x to dataset for future predictions + simulations
      d_cf_mvpa[c("e1_1", "e1_2", "e1_3")] <- d_cf_mvpa_e1[c("e1_1", "e1_2", "e1_3")]
      
      # Predict probability of l3_x (time-varying confounders w3) with assigned e1_x
      d_cf_mvpa$l3_1 <- rbinom(nrow(d_cf_mvpa), 1, predict(m_L3_1, newdata = d_cf_mvpa, type = "response"))
      d_cf_mvpa$l3_2 <- rbinom(nrow(d_cf_mvpa), 1, predict(m_L3_2, newdata = d_cf_mvpa, type = "response"))
      d_cf_mvpa$l3_3 <- rbinom(nrow(d_cf_mvpa), 1, predict(m_L3_3, newdata = d_cf_mvpa, type = "response"))
      d_cf_mvpa$l3_4 <- rbinom(nrow(d_cf_mvpa), 1, predict(m_L3_4, newdata = d_cf_mvpa, type = "response"))
      d_cf_mvpa$l3_5 <- rbinom(nrow(d_cf_mvpa), 1, predict(m_L3_5, newdata = d_cf_mvpa, type = "response"))
      
      p_mat <- predict(m_L3_6, newdata = d_cf_mvpa, type = "probs")
      d_cf_mvpa$l3_6 <- apply(p_mat, 1, function(p) {
        sample(1:4, size = 1, prob = p)})
      d_cf_mvpa$l3_6 <- as.factor(d_cf_mvpa$l3_6)
      
      # Set e3_x - determined from path-specific intervention (distribution of specific MB in a=1 group is equalized to those in a=0 group)
      # Set entire population to unexposed for specific MB (a1=0 [MVPA], a2=observed [screen time], a3=observed [sleep], needed for multivariate exposure model)
      d_cf_mvpa_e3 <- d_cf_mvpa
      d_cf_mvpa_e3$a_1 <- 0
      
      # Predict probabilities of e3_x (exposure W3) and simulate their values
      # Generate error terms for each group's natural outcomes
      eps_cf_mvpa_E3 <- rmvnorm(nrow(d_cf_mvpa_e3), mean = c(0,0,0), sigma = cov_mat_E3)
      d_cf_mvpa_e3$e3_1 <- ifelse(predict(m_E3, newdata = d_cf_mvpa_e3, eq = 1) + eps_cf_mvpa_E3[,1] > 0, 1, 0)
      d_cf_mvpa_e3$e3_2 <- ifelse(predict(m_E3, newdata = d_cf_mvpa_e3, eq = 2) + eps_cf_mvpa_E3[,2] > 0, 1, 0)
      d_cf_mvpa_e3$e3_3 <- ifelse(predict(m_E3, newdata = d_cf_mvpa_e3, eq = 3) + eps_cf_mvpa_E3[,3] > 0, 1, 0)
      
      # Add simulated e3_x to dataset for future predictions + simulations
      d_cf_mvpa[c("e3_1", "e3_2", "e3_3")] <- d_cf_mvpa_e3[c("e3_1", "e3_2", "e3_3")]
      
      # Predict probability of y (outcome) with simulated e1_x and e3_x, and simulated l3_x and simulate y
      d_cf_mvpa$y_cf_mvpa <- rbinom(nrow(d_cf_mvpa), 1, predict(m_Y, newdata = d_cf_mvpa, type = "response"))
      
      # Add simulated cf scenario outcome to svydesign object
      decomp_svydesign_cf_mvpa <- update(
        decomp_svydesign,
        y_cf_mvpa = d_cf_mvpa$y_cf_mvpa
      )
      
      # Compute CDM (RD and RR) for cf scenario = use y_cf
      final_model_cf_mvpa <- svyglm(y_cf_mvpa ~ age_group_wave1 + sex, family = quasibinomial(link = "logit"), design = decomp_svydesign_cf_mvpa)
      risks_cf_mvpa <- svypredmeans(adjustmodel = final_model_cf_mvpa, groupfactor = ~ a)
      RD_cf_mvpa_contrast <- svycontrast(risks_cf_mvpa, quote(`1` - `0`))
      RR_cf_mvpa_contrast <- svycontrast(risks_cf_mvpa, quote(`1` / `0`))
      RD_cf_mvpa_coef_sim_eq[k, i] <- as.numeric(coef(RD_cf_mvpa_contrast)) 
      RR_cf_mvpa_coef_sim_eq[k, i] <- as.numeric(coef(RR_cf_mvpa_contrast))
      
      # Compute group-specific and population-wide risks
      p0_cf_mvpa_sim_eq[k, i] <- as.numeric(coef(risks_cf_mvpa)[1])
      p1_cf_mvpa_sim_eq[k, i] <- as.numeric(coef(risks_cf_mvpa)[2])
      
      pred_probs_cf_mvpa <- predict(final_model_cf_mvpa, type = "response")  
      pop_risk_cf_mvpa <- svymean(pred_probs_cf_mvpa, design = decomp_svydesign_cf_mvpa)
      pop_risk_cf_mvpa_sim_eq[k, i] <- as.numeric(coef(pop_risk_cf_mvpa)[1])
    }
    
    # PATH-SPECIFIC (SLEEP)
    
    for (i in 1:n_mc) {  
      
      # Create dataset for counterfactual/intervention scenario
      d_cf_sleep <- svy_data %>% mutate(across(all_of(vars_to_na), ~ NA))
      
      # Set e1_x - determined from path-specific intervention (distribution of specific MB in a=1 group is equalized to those in a=0 group)
      # Set entire population to unexposed for specific MB (a1=observed [MVPA], a2=0 [sleep], a3=observed [screen time], needed for multivariate exposure model)
      d_cf_sleep_e1 <- d_cf_sleep
      d_cf_sleep_e1$a_2 <- 0
      
      # Predict probability of e1_x (exposure W1) and simulate value
      # Generate error terms for each group's natural/counterfactual outcomes
      eps_cf_sleep_E1 <- rmvnorm(nrow(d_cf_sleep_e1), mean = c(0,0,0), sigma = cov_mat_E1)
      d_cf_sleep_e1$e1_1 <- ifelse(predict(m_E1, newdata = d_cf_sleep_e1, eq = 1) + eps_cf_sleep_E1[,1] > 0, 1, 0)
      d_cf_sleep_e1$e1_2 <- ifelse(predict(m_E1, newdata = d_cf_sleep_e1, eq = 2) + eps_cf_sleep_E1[,2] > 0, 1, 0)
      d_cf_sleep_e1$e1_3 <- ifelse(predict(m_E1, newdata = d_cf_sleep_e1, eq = 3) + eps_cf_sleep_E1[,3] > 0, 1, 0)
      
      # Add simulated e1_x to dataset for future predictions + simulations
      d_cf_sleep[c("e1_1", "e1_2", "e1_3")] <- d_cf_sleep_e1[c("e1_1", "e1_2", "e1_3")]
      
      # Predict probability of l3_x (time-varying confounders w3) with assigned e1_x
      d_cf_sleep$l3_1 <- rbinom(nrow(d_cf_sleep), 1, predict(m_L3_1, newdata = d_cf_sleep, type = "response"))
      d_cf_sleep$l3_2 <- rbinom(nrow(d_cf_sleep), 1, predict(m_L3_2, newdata = d_cf_sleep, type = "response"))
      d_cf_sleep$l3_3 <- rbinom(nrow(d_cf_sleep), 1, predict(m_L3_3, newdata = d_cf_sleep, type = "response"))
      d_cf_sleep$l3_4 <- rbinom(nrow(d_cf_sleep), 1, predict(m_L3_4, newdata = d_cf_sleep, type = "response"))
      d_cf_sleep$l3_5 <- rbinom(nrow(d_cf_sleep), 1, predict(m_L3_5, newdata = d_cf_sleep, type = "response"))
      
      p_mat <- predict(m_L3_6, newdata = d_cf_sleep, type = "probs")
      d_cf_sleep$l3_6 <- apply(p_mat, 1, function(p) {
        sample(1:4, size = 1, prob = p)})
      d_cf_sleep$l3_6 <- as.factor(d_cf_sleep$l3_6)
      
      # Set e3_x - determined from path-specific intervention (distribution of specific MB in a=1 group is equalized to those in a=0 group)
      # Set entire population to unexposed for specific MB (a1=0 [MVPA], a2=observed [screen time], a3=observed [sleep], needed for multivariate exposure model)
      d_cf_sleep_e3 <- d_cf_sleep
      d_cf_sleep_e3$a_2 <- 0
      
      # Predict probabilities of e3_x (exposure W3) and simulate their values
      # Generate error terms for each group's natural/counterfactual outcomes
      eps_cf_sleep_E3 <- rmvnorm(nrow(d_cf_sleep_e3), mean = c(0,0,0), sigma = cov_mat_E3)
      d_cf_sleep_e3$e3_1 <- ifelse(predict(m_E3, newdata = d_cf_sleep_e3, eq = 1) + eps_cf_sleep_E3[,1] > 0, 1, 0)
      d_cf_sleep_e3$e3_2 <- ifelse(predict(m_E3, newdata = d_cf_sleep_e3, eq = 2) + eps_cf_sleep_E3[,2] > 0, 1, 0)
      d_cf_sleep_e3$e3_3 <- ifelse(predict(m_E3, newdata = d_cf_sleep_e3, eq = 3) + eps_cf_sleep_E3[,3] > 0, 1, 0)
      
      # Add simulated e3_x to dataset for future predictions + simulations
      d_cf_sleep[c("e3_1", "e3_2", "e3_3")] <- d_cf_sleep_e3[c("e3_1", "e3_2", "e3_3")]
      
      # Predict probability of y (outcome) with simulated e1_x and e3_x, and simulated l3_x and simulate y
      d_cf_sleep$y_cf_sleep <- rbinom(nrow(d_cf_sleep), 1, predict(m_Y, newdata = d_cf_sleep, type = "response"))
      
      # Add simulated cf scenario outcome to svydesign object
      decomp_svydesign_cf_sleep <- update(
        decomp_svydesign,
        y_cf_sleep = d_cf_sleep$y_cf_sleep
      )
      
      # Compute CDM (RD and RR) for cf scenario = use y_cf
      final_model_cf_sleep <- svyglm(y_cf_sleep ~ age_group_wave1 + sex, family = quasibinomial(link = "logit"), design = decomp_svydesign_cf_sleep)
      risks_cf_sleep <- svypredmeans(adjustmodel = final_model_cf_sleep, groupfactor = ~ a)
      RD_cf_sleep_contrast <- svycontrast(risks_cf_sleep, quote(`1` - `0`))
      RR_cf_sleep_contrast <- svycontrast(risks_cf_sleep, quote(`1` / `0`))
      RD_cf_sleep_coef_sim_eq[k, i] <- as.numeric(coef(RD_cf_sleep_contrast)) 
      RR_cf_sleep_coef_sim_eq[k, i] <- as.numeric(coef(RR_cf_sleep_contrast)) 
      
      # Compute group-specific and population-wide risks
      p0_cf_sleep_sim_eq[k, i] <- as.numeric(coef(risks_cf_sleep)[1])
      p1_cf_sleep_sim_eq[k, i] <- as.numeric(coef(risks_cf_sleep)[2])
      
      pred_probs_cf_sleep <- predict(final_model_cf_sleep, type = "response")  
      pop_risk_cf_sleep <- svymean(pred_probs_cf_sleep, design = decomp_svydesign_cf_sleep)
      pop_risk_cf_sleep_sim_eq[k, i] <- as.numeric(coef(pop_risk_cf_sleep)[1])
    }
    
    # PATH-SPECIFIC (SCREEN TIME)
    
    for (i in 1:n_mc) {  
      
      # Create dataset for counterfactual/intervention scenario
      d_cf_st <- svy_data %>% mutate(across(all_of(vars_to_na), ~ NA))
      
      # Set e1_x - determined from path-specific intervention (distribution of specific MB in a=1 group is equalized to those in a=0 group)
      # Set entire population to unexposed for specific MB (a1=observed [MVPA], a2=observed [sleep], a3=0 [screen time], needed for multivariate exposure model)
      d_cf_st_e1 <- d_cf_st
      d_cf_st_e1$a_3 <- 0
      
      # Predict probability of e1_x (exposure W1) and simulate value
      # Generate error terms for each group's natural/counterfactual outcomes
      eps_cf_st_E1 <- rmvnorm(nrow(d_cf_st_e1), mean = c(0,0,0), sigma = cov_mat_E1)
      d_cf_st_e1$e1_1 <- ifelse(predict(m_E1, newdata = d_cf_st_e1, eq = 1) + eps_cf_st_E1[,1] > 0, 1, 0)
      d_cf_st_e1$e1_2 <- ifelse(predict(m_E1, newdata = d_cf_st_e1, eq = 2) + eps_cf_st_E1[,2] > 0, 1, 0)
      d_cf_st_e1$e1_3 <- ifelse(predict(m_E1, newdata = d_cf_st_e1, eq = 3) + eps_cf_st_E1[,3] > 0, 1, 0)
      
      # Add simulated e1_x to dataset for future predictions + simulations
      d_cf_st[c("e1_1", "e1_2", "e1_3")] <- d_cf_st_e1[c("e1_1", "e1_2", "e1_3")]
      
      # Predict probability of l3_x (time-varying confounders w3) with assigned e1_x
      d_cf_st$l3_1 <- rbinom(nrow(d_cf_st), 1, predict(m_L3_1, newdata = d_cf_st, type = "response"))
      d_cf_st$l3_2 <- rbinom(nrow(d_cf_st), 1, predict(m_L3_2, newdata = d_cf_st, type = "response"))
      d_cf_st$l3_3 <- rbinom(nrow(d_cf_st), 1, predict(m_L3_3, newdata = d_cf_st, type = "response"))
      d_cf_st$l3_4 <- rbinom(nrow(d_cf_st), 1, predict(m_L3_4, newdata = d_cf_st, type = "response"))
      d_cf_st$l3_5 <- rbinom(nrow(d_cf_st), 1, predict(m_L3_5, newdata = d_cf_st, type = "response"))
      
      p_mat <- predict(m_L3_6, newdata = d_cf_st, type = "probs")
      d_cf_st$l3_6 <- apply(p_mat, 1, function(p) {
        sample(1:4, size = 1, prob = p)})
      d_cf_st$l3_6 <- as.factor(d_cf_st$l3_6)
      
      # Set e3_x - determined from path-specific intervention (distribution of specific MB in a=1 group is equalized to those in a=0 group)
      # Set entire population to unexposed for specific MB (a1=observed [MVPA], a2=observed [screen time], a3=0 [sleep], needed for multivariate exposure model)
      d_cf_st_e3 <- d_cf_st
      d_cf_st_e3$a_3 <- 0
      
      # Predict probabilities of e3_x (exposure W3) and simulate their values
      # Generate error terms for each group's natural/counterfactual outcomes
      eps_cf_st_E3 <- rmvnorm(nrow(d_cf_st_e3), mean = c(0,0,0), sigma = cov_mat_E3)
      d_cf_st_e3$e3_1 <- ifelse(predict(m_E3, newdata = d_cf_st_e3, eq = 1) + eps_cf_st_E3[,1] > 0, 1, 0)
      d_cf_st_e3$e3_2 <- ifelse(predict(m_E3, newdata = d_cf_st_e3, eq = 2) + eps_cf_st_E3[,2] > 0, 1, 0)
      d_cf_st_e3$e3_3 <- ifelse(predict(m_E3, newdata = d_cf_st_e3, eq = 3) + eps_cf_st_E3[,3] > 0, 1, 0)
      
      # Add simulated e3_x to dataset for future predictions + simulations
      d_cf_st[c("e3_1", "e3_2", "e3_3")] <- d_cf_st_e3[c("e3_1", "e3_2", "e3_3")]
      
      # Predict probability of y (outcome) with simulated e1_x and e3_x, and simulated l3_x and simulate y
      d_cf_st$y_cf_st <- rbinom(nrow(d_cf_st), 1, predict(m_Y, newdata = d_cf_st, type = "response"))
      
      # Add simulated cf scenario outcome to svydesign object
      decomp_svydesign_cf_st <- update(
        decomp_svydesign,
        y_cf_st = d_cf_st$y_cf_st
      )
      
      # Compute CDM (RD and RR) for cf scenario = use y_cf
      final_model_cf_st <- svyglm(y_cf_st ~ age_group_wave1 + sex, family = quasibinomial(link = "logit"), design = decomp_svydesign_cf_st)
      risks_cf_st <- svypredmeans(adjustmodel = final_model_cf_st, groupfactor = ~ a)
      RD_cf_st_contrast <- svycontrast(risks_cf_st, quote(`1` - `0`))
      RR_cf_st_contrast <- svycontrast(risks_cf_st, quote(`1` / `0`))
      RD_cf_st_coef_sim_eq[k, i] <- as.numeric(coef(RD_cf_st_contrast)) 
      RR_cf_st_coef_sim_eq[k, i] <- as.numeric(coef(RR_cf_st_contrast)) 
      
      # Compute group-specific and population-wide risks
      p0_cf_st_sim_eq[k, i] <- as.numeric(coef(risks_cf_st)[1])
      p1_cf_st_sim_eq[k, i] <- as.numeric(coef(risks_cf_st)[2])
      
      pred_probs_cf_st <- predict(final_model_cf_st, type = "response")  
      pop_risk_cf_st <- svymean(pred_probs_cf_st, design = decomp_svydesign_cf_st)
      pop_risk_cf_st_sim_eq[k, i] <- as.numeric(coef(pop_risk_cf_st)[1])
    }
    
    # MODEL MISSPECIFICATION CHECK: Compare nc to observed data
    model_ob <- svyglm(y ~ age_group_wave1 + sex,
                       family = quasibinomial(link = "logit"),
                       design = decomp_svydesign)
    
    risks_ob <- svypredmeans(adjustmodel = model_ob, groupfactor = ~ a)
    RD_ob_contrast <- svycontrast(risks_ob, quote(`1` - `0`))
    RR_ob_contrast <- svycontrast(risks_ob, quote(`1` / `0`))
    RD_ob_coef_sim[k] <- as.numeric(coef(RD_ob_contrast)) 
    RR_ob_coef_sim[k] <- as.numeric(coef(RR_ob_contrast)) 
    
    # Compute group-specific and population-wide risks
    p0_ob_sim[k] <- as.numeric(coef(risks_ob)[1])
    p1_ob_sim[k] <- as.numeric(coef(risks_ob)[2])
    
    pred_probs_ob <- predict(model_ob, type = "response")  
    pop_risk_ob <- svymean(pred_probs_ob, design = decomp_svydesign)
    pop_risk_ob_sim[k] <- as.numeric(coef(pop_risk_ob)[1])
    
  }
  
  # Collect results
  results <- list(
    # NATURAL COURSE AND OBSERVED
    RD_nc             = RD_nc_coef_sim,
    RD_ob             = RD_ob_coef_sim,
    # ELIMINATION
    RD_cf_joint_el    = RD_cf_joint_coef_sim_el,
    RD_cf_mvpa_el     = RD_cf_mvpa_coef_sim_el,
    RD_cf_sleep_el    = RD_cf_sleep_coef_sim_el,
    RD_cf_st_el       = RD_cf_st_coef_sim_el,
    # EQUALIZATION
    RD_cf_joint_eq    = RD_cf_joint_coef_sim_eq,
    RD_cf_mvpa_eq     = RD_cf_mvpa_coef_sim_eq,
    RD_cf_sleep_eq    = RD_cf_sleep_coef_sim_eq,
    RD_cf_st_eq       = RD_cf_st_coef_sim_eq,
    # NATURAL COURSE AND OBSERVED
    RR_nc             = RR_nc_coef_sim,
    RR_ob             = RR_ob_coef_sim,
    # ELIMINATION
    RR_cf_joint_el    = RR_cf_joint_coef_sim_el,
    RR_cf_mvpa_el     = RR_cf_mvpa_coef_sim_el,
    RR_cf_sleep_el    = RR_cf_sleep_coef_sim_el,
    RR_cf_st_el       = RR_cf_st_coef_sim_el,
    # EQUALIZATION
    RR_cf_joint_eq    = RR_cf_joint_coef_sim_eq,
    RR_cf_mvpa_eq     = RR_cf_mvpa_coef_sim_eq,
    RR_cf_sleep_eq    = RR_cf_sleep_coef_sim_eq,
    RR_cf_st_eq       = RR_cf_st_coef_sim_eq,
    
    # NATURAL COURSE AND OBSERVED
    p0_nc                   = p0_nc_sim,
    p1_nc                   = p1_nc_sim,
    pop_risk_nc             = pop_risk_nc_sim,
    p0_ob                   = p0_ob_sim,
    p1_ob                   = p1_ob_sim,
    pop_risk_ob             = pop_risk_ob_sim,
    # ELIMINATION
    p0_cf_joint_el          = p0_cf_joint_sim_el,
    p1_cf_joint_el          = p1_cf_joint_sim_el ,
    pop_risk_cf_joint_el    = pop_risk_cf_joint_sim_el,
    p0_cf_mvpa_el           = p0_cf_mvpa_sim_el,
    p1_cf_mvpa_el           = p1_cf_mvpa_sim_el,
    pop_risk_cf_mvpa_el     = pop_risk_cf_mvpa_sim_el,
    p0_cf_sleep_el          = p0_cf_sleep_sim_el,
    p1_cf_sleep_el          = p1_cf_sleep_sim_el,
    pop_risk_cf_sleep_el    = pop_risk_cf_sleep_sim_el,
    p0_cf_st_el             = p0_cf_st_sim_el,
    p1_cf_st_el             = p1_cf_st_sim_el,
    pop_risk_cf_st_el       = pop_risk_cf_st_sim_el,
    # EQUALIZATION
    p0_cf_joint_eq          = p0_cf_joint_sim_eq,
    p1_cf_joint_eq          = p1_cf_joint_sim_eq ,
    pop_risk_cf_joint_eq    = pop_risk_cf_joint_sim_eq,
    p0_cf_mvpa_eq           = p0_cf_mvpa_sim_eq,
    p1_cf_mvpa_eq           = p1_cf_mvpa_sim_eq,
    pop_risk_cf_mvpa_eq     = pop_risk_cf_mvpa_sim_eq,
    p0_cf_sleep_eq          = p0_cf_sleep_sim_eq,
    p1_cf_sleep_eq          = p1_cf_sleep_sim_eq,
    pop_risk_cf_sleep_eq    = pop_risk_cf_sleep_sim_eq,
    p0_cf_st_eq             = p0_cf_st_sim_eq,
    p1_cf_st_eq             = p1_cf_st_sim_eq,
    pop_risk_cf_st_eq       = pop_risk_cf_st_sim_eq
  )
  
  # Return results
  return(results)
  
}

# Set parallel computing specifications
library(future)
library(future.apply)

parallel::detectCores()
plan(multisession, workers = 15)

########################################################################################################################################
# RUN FUNCTION
########################################################################################################################################

set.seed(43244)
system.time({
  results_list <- future_lapply(MI_output_final[1:30], function(imp) {
    int_scenarios(data_imp = imp, n_mc = , bs = , n_boot = , sep = '', outcome = '')
  }, future.seed = TRUE)
})