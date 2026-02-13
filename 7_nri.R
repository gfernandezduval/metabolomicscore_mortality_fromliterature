# Title: Net Reclassification Improvement for Score models
# Author: Gonzalo Fernández Duval
# Date: 2025-10-13

# Libraries
library(dplyr)
library(survival)
library(nricens) 
library(ggplot2)
library(ggalluvial)
library(reshape2)
library(gridExtra)

# Data with scores
load("./output_data/df_score_priori.RData")

# New data
df_nri <- df_score_priori

# Covariates
covariates <- c("age10", "sex", "glucose0", "wth_cat", 
                     "smoking0", "alcohol_cat", "hyperten0", 
                     "educ2", "fam_history", "dyslip0", "bmi_cat",
                     "energyt", "getotal", "p14", "group_int")

# Formulas
formula_noscore <- paste("Surv(follow_up19, death19) ~",
                         paste(covariates, collapse = " + "))

formula_decile <- paste("Surv(follow_up19, death19) ~", "score_stand_dec", "+", 
                        paste(covariates, collapse = " + "))

formula_ridge <- paste("Surv(follow_up19, death19) ~", "score_stand_ridge", "+", 
                        paste(covariates, collapse = " + "))


# Models
model_noscore <- coxph(as.formula(formula_noscore), data = df_nri)
model_decile <- coxph(as.formula(formula_decile), data = df_nri)
model_ridge <- coxph(as.formula(formula_ridge), data = df_nri)

# Calculate risk predictions
df_nri$risk_noscore <- predict(model_noscore, type = "risk")
df_nri$risk_decile <- predict(model_decile, type = "risk")
df_nri$risk_ridge <- predict(model_ridge, type = "risk")

#=================================================
#================== TERTILE ======================
#=================================================

# Create tertiles based on the baseline model (noscore)
tertiles <- quantile(df_nri$risk_noscore, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
print("Tertile cutpoints:")
print(tertiles)

# Calculate NRI
nri_tertile_decile <- nricens(time = df_nri$follow_up19,
                              event = df_nri$death19,
                              p.std = df_nri$risk_noscore,
                              p.new = df_nri$risk_decile,
                              t0 = 10,
                              cut = tertiles[-c(1,4)])  # Remove 0 and 1 quantiles


nri_tertile_ridge <- nricens(time = df_nri$follow_up19,
                             event = df_nri$death19,
                             p.std = df_nri$risk_noscore,
                             p.new = df_nri$risk_ridge,
                             t0 = 10,
                             cut = tertiles[-c(1,4)])  

#=================================================
#================== QUINTILE =====================
#=================================================

quintiles <- quantile(df_nri$risk_noscore, probs = c(0, 1/5, 2/5, 3/5, 4/5, 1), na.rm = TRUE)
print("Quintile cutpoints:")
print(quintiles)

# Calculate NRI
nri_quintile_decile <- nricens(time = df_nri$follow_up19,
                             event = df_nri$death19,
                             p.std = df_nri$risk_noscore,
                             p.new = df_nri$risk_decile,
                             t0 = 10,
                             cut = quintiles[-c(1,6)])  # Remove 0 and 1 quantiles


nri_quintile_ridge <- nricens(time = df_nri$follow_up19,
                            event = df_nri$death19,
                            p.std = df_nri$risk_noscore,
                            p.new = df_nri$risk_ridge,
                            t0 = 10,
                            cut = quintiles[-c(1,6)])  

#=======================================================
#================= CONTINUOUS NRI ======================
#=======================================================

# NRI continuous between base model and decile-based model
nri_continuous_decile  <- nricens(time  = df_nri$follow_up19,
                           event = df_nri$death19,
                           p.std = df_nri$risk_noscore,
                           p.new = df_nri$risk_decile,
                           t0    = 10,
                           updown = "diff",
                           cut = 0)       

# NRI continuous between base model and ridge-based model
nri_continuous_ridge <- nricens(time  = df_nri$follow_up19,
                          event = df_nri$death19,
                          p.std = df_nri$risk_noscore,
                          p.new = df_nri$risk_ridge,
                          t0    = 10,
                          updown = "diff",
                          cut   = 0)

