# Title: Replication of PREDIMED scores on PREDIMED at 1 year of follow-up
# Authors: Gonzalo Fernandez-Duval
# Date: 2025-10-13

# Libraries
library(dplyr)
library(survival)
library(survminer)
library(readr)
library(readxl)
library(lubridate)
library(rcompanion)
library(stringr)

#=======================================================
#======================== DATA ========================= 
#=======================================================

# Load needed data
load("./input_data/df_predimed_1year.RData")

# Coef from PREDIMED
predimed_score <- read_xlsx(path = "./input_data/replicate.xlsx") %>%
  select(-n, -availability_nhs_hpfs)

# Order and coeff
order_met <- predimed_score$metabolite

# Names of metabolites
protective_met <- c("Histidine" = "histidine1",
                    "Tryptophan" = "tryptophan1",
                    "SM 24:0" = "c240sm1",      
                    "Leucine" = "leucine1",      
                    "Asparagine" = "asparagine1",      
                    "SM 22:0" = "c220sm1",      
                    "Uridine" = "uridine1",      
                    "Homoarginine" = "homoarginine1" )

risk_met <- c("CAR 4:0;OH" = "c4ohcarnitine1",   
              "Malate" = "malate1",       
              "Cer 18:1/16:0" = "c160ceramided1811",
              "SM 16:0" = "c160sm1",         
              "Dimethylglycine" = "dimethylglycine1", 
              "Cotinine" = "cotinine1",        
              "Beta-hydroxybutyrate" = "betahydroxybutyrate1", 
              "1-Methyladenosine" = "methyladenosine1")

total_met <- c(protective_met, risk_met)

# Number of each metabolite
n_protective <- length(protective_met)
n_risk <- length(risk_met)
n_met <- length(total_met)

# Data for PREDIMED replication
replication_predimed1 <- df_predimed_1year %>%
  select(id,
         all_of(total_met),   
         age0:follow_up19_1) %>% 
  select(id,
         all_of(order_met), 
         age0:follow_up19_1) %>% 
  mutate_at(vars(2:(n_protective + 1)), list(dec = ~11 - ntile(., 10))) %>%       # Inverted deciles (protective)
  mutate_at(vars((n_protective + 2):(n_met + 1)), list(dec = ~ ntile(., 10)))     # Deciles (risk)

#=========================================================
#======================== SCORES ========================= 
#=========================================================

# Score 1
df_score1 <- replication_predimed1 %>%
  mutate(score1 = rowSums(across(ends_with("_dec")), na.rm = T),  
         score1_q = ntile(score1, 5),
         score1_q = factor(score1_q, levels = c("1", "2", "3", "4","5")),
         score1_stand = as.numeric(scale(score1))) %>%
  select(id, score1:score1_stand)
  
# Information for Score 3
score3_coeff <- predimed_score$score_ridge
data_score3 <- replication_predimed1
data_score3[, 2:(n_met + 1)] <- sweep(replication_predimed1[, 2:(n_met + 1)], 2, score3_coeff, "*")  

# Calculate Score 3
df_score3 <- data_score3 %>%
  mutate(score3 = rowSums(across(2:(n_met + 1)), na.rm = T),
         score3_q = ntile(score3, 5),
         score3_q = factor(score3_q, levels = c("1", "2", "3", "4","5")),
         score3_stand = as.numeric(scale(score3))) %>%
  select(id, score3, score3_q, score3_stand)

# Final data replication PREDIMED at 1 year
df_replication_predimed1 <- replication_predimed1 %>%
  left_join(df_score1, by = "id") %>%
  left_join(df_score3, by = "id") 

# Data with scores
save(df_replication_predimed1, file = "./output_data/df_replication_predimed1.RData")

#=====================================================
#================== Cox Models ======================
#===================================================== 

# Define covariate sets
basic_covariates <- c("age10_1", "sex")

full_covariates <- c("age10_1", "sex", "glucose1", "wth_cat_1", 
                     "smoking1", "alcohol_cat1", "hyperten0", 
                     "educ2", "fam_history", "dyslip0", "bmi_cat_1",
                     "energyt1", "getotal1", "p14_1", "group_int")

# Function to create formula string
create_formula <- function(score_var, covariates = NULL) {
  if (is.null(covariates)) {
    paste("Surv(follow_up19_1, death19) ~", score_var)
  } 
  else {
    paste("Surv(follow_up19_1, death19) ~", score_var, "+", 
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

# FINAL MODELS
models <- list(score1_q = create_formula("score1_q"),
               score1_stand = create_formula("score1_stand"),
               score1_q_agesex = create_formula("score1_q", basic_covariates),
               score1_stand_agesex = create_formula("score1_stand", basic_covariates),
               score1_q_full = create_formula("score1_q", full_covariates),
               score1_stand_full = create_formula("score1_stand", full_covariates),
               score3_q = create_formula("score3_q"),
               score3_stand = create_formula("score3_stand"),
               score3_q_agesex = create_formula("score3_q", basic_covariates),
               score3_stand_agesex = create_formula("score3_stand", basic_covariates),
               score3_q_full = create_formula("score3_q", full_covariates),
               score3_stand_full = create_formula("score3_stand", full_covariates))

# Analyze all models
results <- list()
for(model_name in names(models)) {
  results[[model_name]] <- analyze_cox_model(
    formula = models[[model_name]],
    data = df_replication_predimed1,
    model_name = model_name)
}


# MODEL SCORE 1
results$score1_q$summary
results$score1_q_agesex$summary
results$score1_q_full$summary

results$score1_stand$summary
results$score1_stand_agesex$summary
results$score1_stand_full$summary

# MODEL SCORE 3
results$score3_q$summary
results$score3_q_agesex$summary
results$score3_q_full$summary

results$score3_stand$summary
results$score3_stand_agesex$summary
results$score3_stand_full$summary

# To test p for trend
library(car)
linearHypothesis(results$score3_q_full$model, 
                 hypothesis.matrix = c("score1_q2 = 1", 
                                       "score1_q3 = 2", 
                                       "score1_q4 = 3", 
                                       "score1_q5 = 4"))

linearHypothesis(results$score3_q_full$model, 
                 hypothesis.matrix = c("score3_q2 = 1", 
                                       "score3_q3 = 2", 
                                       "score3_q4 = 3", 
                                       "score3_q5 = 4"))
