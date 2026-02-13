# Title: Association between all-cause mortality and different scores alternatives
# Author: Gonzalo Fernández Duval
# Date: 2025-10-13

# Libraries
library(dplyr)
library(survival)
library(survminer)

# Data with scores
load("./output_data/df_score_priori.RData")

#=====================================================
#================== Cox Models ======================
#===================================================== 

# Define covariate sets
basic_covariates <- c("age10", "sex")

full_covariates <- c("age10", "sex", "glucose0", "wth_cat", 
                     "smoking0", "alcohol_cat", "hyperten0", 
                     "educ2", "fam_history", "dyslip0", "bmi_cat",
                     "energyt", "getotal", "p14", "group_int")

# Function to create formula string
create_formula <- function(score_var, covariates = NULL) {
  if (is.null(covariates)) {
    paste("Surv(follow_up19, death19) ~", score_var)
  } 
  else {
    paste("Surv(follow_up19, death19) ~", score_var, "+", 
          paste(covariates, collapse = " + "))
  }
}

# Function to fit and analyze Cox model
analyze_cox_model <- function(formula, data, model_name = "") {
  # Fit model
  model <- coxph(as.formula(formula), data = data)
  
  # Get summary
  model_summary <- summary(model)
  
  # Check proportional hazards assumption
  ph_test <- cox.zph(model)
  
  return(list(model = model,
              summary = model_summary,
              ph_test = ph_test))
}
##### FINAL MODELS ######
models <- list(score_crude1 = create_formula("score_q_dec"),
                     score_agesex1 = create_formula("score_q_dec", basic_covariates),
                     score_quintiles1 = create_formula("score_q_dec", full_covariates),
                     score_crude1_stand = create_formula("score_stand_dec"),
                     score_agesex1_stand = create_formula("score_stand_dec", basic_covariates),
                     score_quintiles1_stand = create_formula("score_stand_dec", full_covariates),
                     score_crude2 = create_formula("score_q_ridge"),
                     score_agesex2 = create_formula("score_q_ridge", basic_covariates),
                     score_quintiles2 = create_formula("score_q_ridge", full_covariates),
                     score_crude2_stand = create_formula("score_stand_ridge"),
                     score_agesex2_stand = create_formula("score_stand_ridge", basic_covariates),
                     score_quintiles2_stand = create_formula("score_stand_ridge", full_covariates)
      )

# Analyze all models
results <- list()
for(model_name in names(models)) {
  results[[model_name]] <- analyze_cox_model(
    formula = models[[model_name]],
    data = df_score_priori,
    model_name = model_name)
}

# MODEL SCORE DECILE
results$score_crude1$summary
results$score_agesex1$summary
results$score_quintiles1$summary

results$score_crude1_stand$summary
results$score_agesex1_stand$summary
results$score_quintiles1_stand$summary

# MODEL SCORE RIDGE
results$score_crude2$summary
results$score_agesex2$summary
results$score_quintiles2$summary

results$score_crude2_stand$summary
results$score_agesex2_stand$summary
results$score_quintiles2_stand$summary


# To test p for trend
library(car)
linearHypothesis(results$score_quintiles1$model, 
                 hypothesis.matrix = c("score_q_dec2 = 1", 
                                       "score_q_dec3 = 2", 
                                       "score_q_dec4 = 3", 
                                       "score_q_dec5 = 4"))

linearHypothesis(results$score_quintiles2$model, 
                 hypothesis.matrix = c("score_q_ridge2 = 1", 
                                       "score_q_ridge3 = 2", 
                                       "score_q_ridge4 = 3", 
                                       "score_q_ridge5 = 4"))

