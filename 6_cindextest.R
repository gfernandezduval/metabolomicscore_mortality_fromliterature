# Title: C-index comparison for baseline PREDIMED
# Author: Gonzalo Fernández Duval
# Date: 2025-10-13

# Library
library(dplyr)
library(survival)
library(compareC)
library(survC)

#===================================
#=========== DATA ================
#===================================

# Load data
load("./output_data/df_score_priori.RData")

# Score
#score <- "score_stand_dec"
score <- "score_stand_ridge"    

# Define covariates
var_modelvar <- c("age10", "sex", "glucose0", "wth_cat", 
                     "smoking0", "alcohol_cat", "hyperten0", 
                     "educ2", "fam_history", "dyslip0", "bmi_cat",
                     "energyt", "getotal", "p14", "group_int")

var_modelscore <- c(score, var_modelvar)

#===================================
#=========== MODELS ================
#===================================

# Formula
form_var <- paste("Surv(follow_up19, death19) ~", paste(var_modelvar, collapse =  " + "))
form_score <- paste("Surv(follow_up19, death19) ~", paste(var_modelscore, collapse =  " + "))

# Models
model_var <- coxph(as.formula(form_var), data = df_score_priori)
summary(model_var)
model_score <- coxph(as.formula(form_score), data = df_score_priori)
summary(model_score)

# C-index with IC
c_var <- cindex_calc(model_var)
c_score <- cindex_calc(model_score)

# Prediction
pred_var <- predict(model_var, type = "lp")
pred_score <- predict(model_score, type = "lp")

# C-index comparison test
cindex_test <- compareC(timeX = df_score_priori$follow_up19,
                        statusX = df_score_priori$death19,
                        scoreY = pred_score,
                        scoreZ = pred_var)

cindex_test
