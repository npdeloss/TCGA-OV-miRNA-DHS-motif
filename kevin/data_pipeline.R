# import libraries
require(data.table)
require(dtplyr)
require(tidyverse)
require(randomForestSRC)
require(survival)

# user defined libraries
source("TCGA PFS.R")

# read in drug data ####
drug_df <- read.table("data/clinical.patient.drug.txt", header=T)

# read in clinical df
clinical_df <- read.table("data/clinical.patient.txt", header=T)

# merge in PFS
clinical_df <- merge(clinical_df, TCGA_PFS(drug_df, clinical_df), all.x=TRUE, all.y=FALSE)


# recode stage
clinical_df$stage_recode[clinical_df$stage_event_clinical_stage 
                         %in% c("Stage IA", "Stage IB", "Stage IC")] <- 1
clinical_df$stage_recode[clinical_df$stage_event_clinical_stage 
                         %in% c("Stage IIA", "Stage IIB", "Stage IIC")] <- 2
clinical_df$stage_recode[clinical_df$stage_event_clinical_stage 
                         %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC")] <- 3
clinical_df$stage_recode[clinical_df$stage_event_clinical_stage 
                         %in% c("IV")] <- 4

clinical_df$stage_recode2[clinical_df$stage_recode == 1] <- "lower stage"
clinical_df$stage_recode2[clinical_df$stage_recode == 2] <- "lower stage"
clinical_df$stage_recode2[clinical_df$stage_recode == 3] <- "higher stage"
clinical_df$stage_recode2[clinical_df$stage_recode == 4] <- "higher stage"

clinical_df$stage_recode2 <- as.factor(clinical_df$stage_recode2)


#  recode OS
# if dead, use days_to_death
# if alive, use days_to_last_followup
dead_idx = clinical_df$vital_status == "Dead"
alive_idx = clinical_df$vital_status == "Alive"

clinical_df$os_time[dead_idx] = clinical_df[dead_idx, "days_to_death"]
clinical_df$os_time[alive_idx] = clinical_df[alive_idx, "days_to_last_followup"]

clinical_df$os_code[clinical_df$vital_status == "Alive"] <- 0
clinical_df$os_code[clinical_df$vital_status == "Dead"] <- 1

# drug tab read in and recode ####
# read in drug table data

# miR read in and recode ####
# read in miR data
miR_df <- fread("data/DHS_subset.miRNA_quantifications.reads_per_million_miRNA_mapped.transposed.txt", head=T)

# average duplicated data
ave_miR_df <- miR_df
ave_miR_df[, miRNA:=NULL]

keys <- names(ave_miR_df)[!grepl('hsa', names(ave_miR_df))]

ave_miR_df <- ave_miR_df[, lapply(.SD, mean), keys]

# merge with clinical file and prepare data
mrg_df <- merge(clinical_df, ave_miR_df)
miRs <- names(mrg_df)[grepl('hsa', names(mrg_df))]
summary(mrg_df$days_to_last_followup)

# impute and regress data ####
trunk_df <- mrg_df[c(miRs, "os_time", "os_code")]
imputed_df <- impute.rfsrc(Surv(os_time, os_code)~., data = trunk_df)
os_rfsrc <- rfsrc(Surv(os_time, os_code)~., data = imputed_df, importance = "permute", seed = 1223)
write.csv(as.data.frame(-sort(os_rfsrc$importance)), "out/os_gene_importance.csv")
summary(survfit(Surv(trunk_df$os_time, trunk_df$os_code)~1), times = 60)
survfit(Surv(trunk_df$os_time, trunk_df$os_code)~1)

# impute and regress PFS ####
trunk_df <- mrg_df[c(miRs, "PFS_days", "PFS_recurr")]
imputed_df <- impute.rfsrc(Surv(PFS_days, PFS_recurr)~., data = trunk_df)
pfs_rfsrc <- rfsrc(Surv(PFS_days, PFS_recurr)~., data = imputed_df, importance = "permute", seed = 1223)
write.csv(as.data.frame(-sort(pfs_rfsrc$importance)), "out/pfs_gene_importance.csv")
summary(survfit(Surv(trunk_df$PFS_days, trunk_df$PFS_recurr)~1), times = 60)
survfit(Surv(trunk_df$PFS_days, trunk_df$PFS_recurr)~1)

# impute and classify stage of presentation ####
trunk_df <- mrg_df[c(miRs, "stage_recode2")]
imputed_df <- impute.rfsrc(stage_recode2~., data = trunk_df)
stage_rfsrc <- rfsrc(stage_recode2~., data = imputed_df, importance = "permute", seed = 1223)
write.csv(as.data.frame(-sort(stage_rfsrc$importance)), "out/stage_gene_importance.csv")
prop.table(trunk_df$stage_recode2)

save.image("miR OvCa.RData")