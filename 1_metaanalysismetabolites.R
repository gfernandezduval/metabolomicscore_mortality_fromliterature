# Title: Meta-analysis metabolite and scores across studies
# Author: Gonzalo Fernández Duval
# Date: 2025-10-13

# Libraries
library(readxl)
library(writexl)
library(dplyr)
library(meta)
library(showtext)
library(dplyr)
library(ggplot2)
library(stringr)

#==========================================================
#========================== DATA ==========================
#==========================================================

############### METABOLITES #####################

# Data from metabolites
metabolites_hr <- read_excel("./output_data/individual_associations_allstudies.xlsx") 

# Metabolites direction
metabolites_direction <- read_excel("./input_data/replicate.xlsx") %>%
  mutate(risk = ifelse(score_decile == "Decile", "risk", "protective")) %>%
  select(metabolite, risk)

# Creating two datasets: one for q5 vs. q1 and the other for standardized values
metabolites_quint <- metabolites_hr %>%
  mutate(comparison = "Q5 vs. Q1") %>%
  select(study, metabolite, comparison, model, hr_q5, lower_ci_q5, upper_ci_q5) %>%
  rename(hr = hr_q5,
         lower_ci = lower_ci_q5,
         upper_ci = upper_ci_q5) %>%
  mutate(hr_ci = paste0(round(hr, digits = 2), " (", round(lower_ci, digits = 2), ", ", round(upper_ci, digits = 2), ")"))
  

metabolites_stand <- metabolites_hr %>%
  mutate(comparison = "Standardized") %>%
  select(study, metabolite, comparison, model, hr_std, lower_ci_std, upper_ci_std) %>%
  rename(hr = hr_std,
         lower_ci = lower_ci_std,
         upper_ci = upper_ci_std) %>%
  mutate(hr_ci = paste0(round(hr, digits = 2), " (", round(lower_ci, digits = 2), ", ", round(upper_ci, digits = 2), ")"))

########################### FINAL DATASETS ###################################

# Q5 vs. Q1
df_quintiles <- metabolites_quint

# Standardized
df_standardized <- metabolites_stand

#===================================================================
#========================== Meta-Analysis ==========================
#===================================================================

# Function to execute meta-analysis with fixed effects
metaanalysis <- function(data) {
  # Convert HR to log(HR) for analysis
  data$log_hr <- log(data$hr)
  
  # Calculate SE of log(HR)
  data$se_log_hr <- (log(data$upper_ci) - log(data$lower_ci)) / (2 * qnorm(0.975))
  
  # Fixed effects meta-analysis
  ma <- metagen(TE = log_hr, 
                seTE = se_log_hr,
                data = data,
                studlab = study,
                fixed = TRUE,
                sm = "HR",
                comb.fixed = TRUE,
                comb.random = FALSE,
                method.tau = "DL")
  
  # Extract results
  hr_fixed <- exp(ma$TE.fixed)
  lower_fixed <- exp(ma$lower.fixed)
  upper_fixed <- exp(ma$upper.fixed)
  hr_ci_fixed <- paste0(round(hr_fixed, 2), " (", 
                        round(lower_fixed, 2), ", ", 
                        round(upper_fixed, 2), ")")
  
  i2 <- round(ma$I2, 3)
  p_het <- round(ma$pval.Q, 3)
  
  return(list(hr_fixed_effect = hr_ci_fixed, 
              i_2 = i2, 
              p = p_het))
}

# Create results table
results_quintile <- data.frame()
results_standardized <- data.frame()

# Obtain all unique combinations of metabolite and model
combinations_quintile <- unique(df_quintiles[, c("metabolite", "model")])
combinations_standardized <- unique(df_standardized[, c("metabolite", "model")])

##################### QUINTILE #########################
# Iterate for each metabolite and model combination
for (i in 1:nrow(combinations_quintile)) {
  metabolite <- combinations_quintile$metabolite[i]
  model <- combinations_quintile$model[i]
  
  # Filter data for the specific combination
  data_subset <- df_quintiles[df_quintiles$metabolite == metabolite & 
                                df_quintiles$model == model, ]
  
  # Extract HR for each study
  hr_predimed <- data_subset$hr_ci[data_subset$study == "PREDIMED"]
  hr_clsa <- data_subset$hr_ci[data_subset$study == "CLSA"]  
  hr_nhshpfs <- data_subset$hr_ci[data_subset$study == "NHS_HPFS"]
  
  # Fill NA's if there is no data available on specific study
  if (length(hr_predimed) == 0) 
    hr_predimed <- NA
  if (length(hr_clsa) == 0) 
    hr_clsa <- NA
  if (length(hr_nhshpfs) == 0) 
    hr_nhshpfs <- NA
  
  # Execute meta-analysis for those who had at least two studies
  if (nrow(data_subset) >= 2) {
    meta_result <- metaanalysis(data_subset)
    hr_fixed_effect <- meta_result$hr_fixed_effect
    i_2 <- meta_result$i_2
    p_value <- meta_result$p}
  else {
    hr_fixed_effect <- NA
    i_2 <- NA
    p_value <- NA}
  
  # Create results row
  row_result <- data.frame(metabolite = metabolite,
                            model = model,
                            hr_predimed = hr_predimed,
                            hr_clsa = hr_clsa,
                            hr_nhshpfs = hr_nhshpfs,
                            hr_fixed_effect = hr_fixed_effect,
                            i_2 = i_2,
                            p = p_value)
  
  # Add to final table
  results_quintile <- rbind(results_quintile, row_result)
}

##################### STANDARDIZED #########################
# Iterate for each metabolite and model combination
for (j in 1:nrow(combinations_standardized)) {
  metabolite_std <- combinations_standardized$metabolite[j]
  model_std <- combinations_standardized$model[j]
  
  # Filter data for the specific combination
  data_subset_std <- df_standardized[df_standardized$metabolite == metabolite_std & 
                                   df_standardized$model == model_std, ]
  
  # Extract HR for each study
  hr_predimed_std <- data_subset_std$hr_ci[data_subset_std$study == "PREDIMED"]
  hr_clsa_std <- data_subset_std$hr_ci[data_subset_std$study == "CLSA"]  
  hr_nhshpfs_std <- data_subset_std$hr_ci[data_subset_std$study == "NHS_HPFS"]
  
  # Fill NA's if there is no data available on specific study
  if (length(hr_predimed_std) == 0) 
    hr_predimed_std <- NA
  if (length(hr_clsa_std) == 0) 
    hr_clsa_std <- NA
  if (length(hr_nhshpfs_std) == 0) 
    hr_nhshpfs_std <- NA
  
  # Execute meta-analysis for those who had at least two studies
  if (nrow(data_subset_std) >= 2) {
    meta_result_std <- metaanalysis(data_subset_std)
    hr_fixed_effect_std <- meta_result_std$hr_fixed_effect
    i_2_std <- meta_result_std$i_2
    p_value_std <- meta_result_std$p}
  else {
    hr_fixed_effect_std <- NA
    i_2_std <- NA
    p_value_std <- NA}
  
  # Crear fila para el resultado
  row_result_std <- data.frame(metabolite = metabolite_std,
                           model = model_std,
                           hr_predimed_std = hr_predimed_std,
                           hr_clsa_std = hr_clsa_std,
                           hr_nhshpfs_std = hr_nhshpfs_std,
                           hr_fixed_effect_std = hr_fixed_effect_std,
                           i_2_std = i_2_std,
                           p_std = p_value_std)
  
  # Agregar a la tabla final
  results_standardized <- rbind(results_standardized, row_result_std)
}

#==========================================================================
#========================== ALL Results and save ==========================
#==========================================================================

# Merge ally
df_all <- results_quintile %>%
  inner_join(results_standardized, by = c("metabolite", "model"))

# Save excel
write_xlsx(df_all, path = "./output_data/meta_analysis_metabolites.xlsx")

# Assign risk
df_all <- df_all %>%
  left_join(metabolites_direction, by = c("metabolite"))

#==========================================================
#========================== PLOT ==========================
#==========================================================

# Data for plot
data_plot <- df_all %>%
  filter(metabolite != "Decile-based score" & metabolite != "Ridge-based score") %>%
  filter(model == "Age and Sex adjusted")

# Font for plot
font_add(family = "ArialNova", 
         regular = "C:/Users/gfern/AppData/Local/Microsoft/Windows/Fonts/ArialNova.ttf")
showtext_auto()
myfont <- "ArialNova"

############# Quintiles
# Extract HR and CI's
df_forest_quint <- data_plot %>%
  mutate(hr_num = as.numeric(str_extract(hr_fixed_effect, "^[0-9\\.]+")),
          ci_low = as.numeric(str_extract(hr_fixed_effect, "(?<=\\().+?(?=,)")),
          ci_up  = as.numeric(str_extract(hr_fixed_effect, "(?<=, ).+?(?=\\))")),
         color_group = ifelse(risk == "risk", "#d25c4d", "#456e9f")) %>%
  arrange(hr_num) %>% 
  mutate(metabolite = factor(metabolite, levels = metabolite))

# Forest plot
forest_quint <- ggplot(df_forest_quint, aes(x = hr_num, y = metabolite)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_up, color = color_group), 
                 height = 0.25, linewidth = 0.9) +
  geom_point(aes(color = color_group), size = 4)  + #shape = 16
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(trans = "log10",
                     breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2),
                     labels = c("0.5", "0.75", "1", "1.25", "1.5", "2")) +
  theme_classic(base_family = myfont) +
  theme(axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = "Hazard Ratio (95% CI) for Q5 vs. Q1 values")

forest_quint

############# STD
# Extract HR and CI's
df_forest_std <- data_plot %>%
  mutate(hr_num_std = as.numeric(str_extract(hr_fixed_effect_std, "^[0-9\\.]+")),
         ci_low = as.numeric(str_extract(hr_fixed_effect_std, "(?<=\\().+?(?=,)")),
         ci_up  = as.numeric(str_extract(hr_fixed_effect_std, "(?<=, ).+?(?=\\))")),
         color_group = ifelse(risk == "risk", "#d25c4d", "#456e9f")) %>%
  arrange(hr_num_std) %>% 
  mutate(metabolite = factor(metabolite, levels = metabolite))

# Forest plot
forest_std <- ggplot(df_forest_std, aes(x = hr_num_std, y = metabolite)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_up, color = color_group), 
                 height = 0.25, linewidth = 0.9) +
  geom_point(aes(color = color_group), size = 4)  + #shape = 16
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(trans = "log10",
                     breaks = c(0.8, 0.9, 1, 1.1, 1.2),
                     labels = c(0.8, 0.9, 1, 1.1, 1.2)) +
  theme_classic(base_family = myfont) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = "Hazard Ratio (95% CI) for +1 SD values")

forest_std

# Save as high-quality PDF
ggsave("./figures/fixedmeta_quintiles.pdf", 
       plot = forest_quint,
       width = 5.5, height = 6, 
       units = "in",
       bg = "white")

# Save as high-quality PDF
ggsave("./figures/fixedmeta_standardized.pdf", 
       plot = forest_std,
       width = 5.5, height = 6, 
       units = "in",
       bg = "white")


