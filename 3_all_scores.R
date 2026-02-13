# Title: All previous Score in one database (Decile-based and Ridge-based)
# Author: Gonzalo Fernández Duval
# Date: 2025-10-13

# Libraries
library(dplyr)
library(haven)

##### All Data
load("./output_data/df_mortality_priori.RData")

# Score using deciles
score_dec <- df_mortality_priori %>% 
  rowwise() %>%                            
  mutate(score_dec = sum(c_across(ends_with("_dec")))) %>%
  ungroup() %>%
  mutate(score_q_dec = as.factor(ntile(score_dec, 5)),
         score_stand_dec = as.numeric(scale(score_dec))) %>%
  as.data.frame() %>%
  select(id, score_dec, score_q_dec, score_stand_dec)

##### Score Ridge
load("./output_data/score_ridge.RData")
score_ridge <- score_ridge %>%
  select(id, score_ridge, score_q_ridge, score_stand_ridge)

#=====================================
#=========== FINAL DATA ==============
#=====================================

# Merge all scores
df_score_priori <- df_mortality_priori %>%
  left_join(score_dec, by = "id") %>%
  left_join(score_ridge, by = "id")

# Save
save(df_score_priori, file = "./output_data/df_score_priori.Rdata")


