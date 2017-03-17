# 3-12-2017
# Kevin Blansit, M.S.
# UC San Diego Biomedical Informatics

# load libraries ####
require(dplyr)

# function ####
TCGA_PFS <- function(drug_df, clinical_df) {
  # INPUTS:
  #  drug_df: the drug biotab
  #  clinical_df: the clinical biotab (needed for last followup)
  
  # primary time ####
  # find end date of adjuvant therapy for mx_therapy_time
  initial <- drug_df %>%
    group_by(bcr_patient_barcode) %>%
    filter(regimen_indication == "ADJUVANT") %>%
    filter(days_to_drug_therapy_start == max(days_to_drug_therapy_start)) %>%
    select(bcr_patient_barcode, mx_therapy_time = days_to_drug_therapy_start)
  
  # find start date for either progression or recurrence therapy for min_recurr_time
  recurr <- drug_df %>%
    group_by(bcr_patient_barcode) %>%
    filter(regimen_indication == "PROGRESSION" | regimen_indication == "RECURRENCE") %>%
    filter(days_to_drug_therapy_end == min(days_to_drug_therapy_end)) %>%
    select(bcr_patient_barcode, min_recurr_time = days_to_drug_therapy_end)
  
  # merge dfs
  mrg_df1 <- merge(initial, recurr, all.x=TRUE, all.y=FALSE)
  
  # determine who has an observed recurrance
  mrg_df1$PFS_recurr <- as.integer(is.na(mrg_df1$min_recurr_time))
  
  # merge clinical
  mrg_df2 <- merge(mrg_df1, clinical_df[c("bcr_patient_barcode", "days_to_last_followup")], all=FALSE)
  
  # determine recurr time for recurr observations
  recurred_idx <- as.logical(!mrg_df2$PFS_recurr)
  mrg_df2$PFS_days[recurred_idx] <- (mrg_df2$min_recurr_time - mrg_df2$mx_therapy_time)[recurred_idx]
  mrg_df2$PFS_days[!recurred_idx] <- mrg_df2$days_to_last_followup[!recurred_idx]
  
  # select and return 
  return(mrg_df2[c("bcr_patient_barcode", "PFS_recurr", "PFS_days")])
  
}