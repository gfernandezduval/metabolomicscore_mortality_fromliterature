# A literature-informed plasma metabolomic signature predicts mortality across several independent cohorts with reproducible performance

The R-code Scripts in this repository cover the main analysis:

-  *0_data.R* : After handling metabolomic data from PREDIMED, indicate the metabolite selection and their direction.
-  *1_metaanalysismetabolites.R* : Meta-analysis of the metabolite selection across studies
-  *2_score_ridge.R* : A priori score, ridge-based version
-  *3_all_scores.R* : Calculate decile-based score and merge the data with ridge-based score.
-  *4_score_mortality.R* : Association of baseline PREDIMED scores with all-cause mortality
-  *5_predimed_1year.R* : Validation of scores in PREDIMED at 1 year of follow-up (example compatible with other replications)
-  *6_cindextest.R* : C-index comparison for baseline PREDIMED (same as PREDIMED at 1 year of follow-up)
-  *7_nri.R* : Continuous Net Reclassification Improvement for both scores

If you have any questions, please contact Gonzalo Fernández-Duval ( ghfernandezd@unav.es )
