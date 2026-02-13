# Title: Data for Score "a priori" associated with all-cause mortality
# Author: Gonzalo Fernández Duval
# Date: 2025-10-13

# Libraries
library(dplyr)
library(readxl)
library(haven)
library(writexl)
library(readr)
library(lubridate)
library(rcompanion)

#===============================
#=========== DATA ==============
#===============================

# We load the worked data with metabolites (already worked)
load(file = ".\\input_data\\df_mortality.Rdata")  

# Formal metabolites names
namesformal <- read_excel(".\\input_data\\metabolitenames.xlsx") 

#=====================================
#=========== PROTECTIVE ==============
#=====================================

# A priori metabolites selected with protective association with all-cause mortality
protective_metabo <- c("Histidine",    # NO in 38 
                       "Tryptophan",   # NO in 38
                       "SM 24:0",      # NO in 38
                       "Leucine",      # NO in 38
                       "Asparagine",   # NO in 38    
                       "SM 22:0",      # NO in 38
                       "Uridine",      # NO in 38 
                       "Homoarginine" # YES in 38   
                       )

# Column names
protective_names <- namesformal %>%
  subset(name %in% protective_metabo) %>%
  select(metabolites) %>% 
  pull(metabolites) 

# Protective metabolites in dataset
protective <- df_mortality %>%
  select(id, 
         all_of(protective_names)) 
  
# Number of variables
n_protective = ncol(protective) - 1

# Final protective data
df_protective <- protective %>%
  mutate_at(vars(2:(n_protective + 1)), list(dec = ~11 - ntile(., 10)))

#===============================
#=========== RISK ==============
#===============================

# A priori metabolites selected with risk association with all-cause mortality
risk_metabo <- c("CAR 4:0;OH",   # YES in 38 
                 "Malate",       # NO in 38
                 "Cer 18:1/16:0",   # NO in 38
                 "SM 16:0",         # NO in 38
                 "Dimethylglycine", # NO in 38
                 "Cotinine",        # Yes in 38
                 "Beta-hydroxybutyrate", # NO in 38
                 "1-Methyladenosine"
                 )

# Column names
risk_names <- namesformal %>%
  subset(name %in% risk_metabo) %>%
  select(metabolites) %>% 
  pull(metabolites) 

# Risk metabolites in dataset
risk <- df_mortality %>%
  select(id, 
         all_of(risk_names)) 

# Number of variables
n_risk = ncol(risk) - 1

# Final risk data
df_risk <- risk %>%
  mutate_at(vars(2:(n_risk + 1)), list(dec = ~ ntile(., 10))) 

#=====================================
#=========== FINAL DATA ==============
#=====================================

# Covariables for multivariable model
adj_var <- c("age10",
             "sex",
             "diabetes0",
             "glucose0",
             "wth_cat",
             "smoking0", 
             "alcohol_cat", 
             "hyperten0", 
             "educ2",
             "fam_history", 
             "dyslip0", 
             "bmi_cat",
             "energyt", 
             "getotal",
             "p14",
             "group_int")

# Adjustment and main variables
adjust_var <- df_mortality %>%
  select(id, 
         all_of(adj_var),
         follow_up19, death19) 

# Merge all the data and calculate final score (Q5 and standardized)
df_mortality_priori <- df_protective %>%
  inner_join(df_risk, by = "id") %>%
  select(-matches("_dec"), everything(), matches("_dec")) %>%
  inner_join(adjust_var, by = "id") 

# Checking NA's
sapply(df_mortality_priori, function(x) sum(is.na(x)))

#Saving information for future Scripts
save(df_mortality_priori, file = ".\\output_data\\df_mortality_priori.Rdata")



