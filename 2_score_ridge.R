# Title: Score a priori association with mortality, using Ridge (all metabolites included)
# Author: Gonzalo Fernández Duval
# Date: 2025-10-13

# Libraries
library(dplyr)
library(survival)
library(glmnet)
library(cvTools)
library(writexl)

# Data
load("./output_data/df_mortality_priori.RData")
names(df_mortality_priori)

# Define metabolite column range
first_m <- as.integer(2)
last_m <- as.integer(17)  

# Prepare modeling data
df_metabolomic <- df_mortality_priori %>%
  select(id,
         all_of(first_m:last_m),       
         follow_up19, 
         death19)

# Create matrix for modeling
x_metabolites <- as.matrix(df_metabolomic[,first_m:last_m])
metabolites_names <- colnames(x_metabolites)

# Data for coefficients to each participant
df_coefficients <- select(df_metabolomic, -follow_up19, -death19)
df_coefficients[,first_m:last_m] <- NA     

#Data to add the metabolite coefficients in each fold
fold_metabolites <- data.frame(metabolite = metabolites_names,
                               matrix(NA, 
                                      nrow = length(metabolites_names), 
                                      ncol = 10, 
                                      dimnames = list(NULL, paste0("fold", 1:10))))

#======================================================
#================= 10-iteration LOFO  =================
#======================================================

# Set seed
set.seed(7102025) 

# Create 10 folds
k = 10
n = nrow(df_metabolomic)
folds = cvFolds(n, K = k)

# Iterate through folds
for (i in 1:k) {
  
  cat("Iteration Fold ", i, '\n')
  
  # Train and test indices
  train_idx <- folds$subsets[folds$which != i]
  test_idx <- folds$subsets[folds$which == i]
  
  # Finding the best lambda
  cv_fit <- cv.glmnet(x = x_metabolites[train_idx, ], 
                      y = Surv(df_metabolomic$follow_up19[train_idx], df_metabolomic$death19[train_idx]), 
                      alpha = 0,   # RIDGE
                      nfolds = 10,
                      family = "cox") 
  
  # Get coefficients
  model_coef <- coef(cv_fit, s= cv_fit$lambda.1se)
  selected_idx <- model_coef@i + 1
  coef_values <- model_coef@x      
  
  # Store selected metabolites and their coefficients for this fold
  fold_metabolites[selected_idx, paste0("fold", i)] <- coef_values
  
  # Save coefficients for test set participants
  test_ids <- df_metabolomic$id[test_idx]
  for(test_id in test_ids) {
    id_row <- which(df_coefficients$id == test_id)
    df_coefficients[id_row, selected_idx + 1] <- coef_values
  }
  
}

#=========================================================
#================= RIDGE a priori Score  =================
#=========================================================

# Data frame with individual multiplications for each metabolite
df_multiplied <- df_coefficients
df_multiplied[, first_m:last_m] <- df_coefficients[, first_m:last_m] * x_metabolites

# Final baseline Score
score_ridge <- df_multiplied %>%
  mutate(score_ridge = rowSums(select(., all_of(first_m):all_of(last_m)), na.rm = T),    # Apply coeff to each metabolite
         score_q_ridge = ntile(score_ridge, 5),
         score_q_ridge = factor(score_q_ridge, levels = c("1", "2", "3", "4", "5")),
         score_stand_ridge = as.numeric(scale(score_ridge)),
         score_ridge = round(score_ridge, digits = 4))

# Final table of each fold 
fold_metabolites_ridge <- fold_metabolites %>%
  mutate(coeff_mean = round(rowMeans(across(fold1:fold10), na.rm = T), digits = 4)) %>%
  arrange(desc(coeff_mean))

# Save results
save(df_coefficients, file = "./output_data/df_coefficients.RData")
