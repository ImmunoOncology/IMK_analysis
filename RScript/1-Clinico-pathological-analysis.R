###################################################
## 
## Clinico pathological analysis
##
## - Descriptive analysis
## - Bivariate analysis 
## - Survival analysis
##
## Manuscript --> Table S1, Table S2, Table S12
##
###################################################

library(survival)
library(survminer)

source("help_fn.R")

clinical <- read.delim("../data/clinical.txt")
clinical$Response <- factor(clinical$Response, levels = c("Bad", "Good"))

vars <- c(
  "Age","Sex","Stage","Diabetes","DL","Hypertension","BRAF","M_Lung","M_node","CNS_metastasis","num_M1","LDH_m1","LDH_treatment_previous_to_IT","Lymphocytes_treatment_previous_to_IT","Neutrophiles_treatment_previous_to_IT","Platelets_treatment_previous_to_IT","No_Previous_lines","Toxicity_IT","Maximum_toxicity_grade"
)

# Melanoma ----------------------------------------------------------------

melanoma_test <- do_bivariate_clinical(clinical)
melanoma_descriptive <- do_descriptive(clinical, vars = melanoma_test$Variable)
melanoma_survival_OS <- do_survival(clinical, vars, time = "OS", event = "OS_event")
melanoma_survival_PFS <- do_survival(clinical, vars, time = "PFS", event = "PFS_event")

# Melanoma cutaneous ------------------------------------------------------

clinica_cutaneous <- clinical[clinical$Diagnosis%in%"cutaneous", ]
cutaneous_test <- do_bivariate_clinical(clinica_cutaneous)
cutaneous_descriptive <- do_descriptive(clinica_cutaneous, vars = melanoma_test$Variable)
cutaneous_survival_OS <- do_survival(clinica_cutaneous, vars, time = "OS", event = "OS_event")
cutaneous_survival_PFS <- do_survival(clinica_cutaneous, vars, time = "PFS", event = "PFS_event")


Table.S1 <- data.frame(
  Variable = melanoma_test$Variable, 
  LRT.melanoma = melanoma_test$p.value, 
  LRT.cutaneous = cutaneous_test$p.value, 
  Logrank.OS.melanoma = melanoma_survival_OS$Logrank[match(melanoma_test$Variable, names(melanoma_survival_OS$variable))], 
  Logrank.OS.cutaneous = cutaneous_survival_OS$Logrank[match(melanoma_test$Variable, names(cutaneous_survival_OS$variable))], 
  Logrank.PFS.melanoma = melanoma_survival_PFS$Logrank[match(melanoma_test$Variable, names(melanoma_survival_PFS$variable))], 
  Logrank.PFS.cutaneous = cutaneous_survival_PFS$Logrank[match(melanoma_test$Variable, names(cutaneous_survival_PFS$variable))]
)

Table.S2 <- data.frame(
  HR.OS.melanoma = melanoma_survival_OS[, -ncol(melanoma_survival_OS)], 
  HR.OS.cutaneous = cutaneous_survival_OS[, -ncol(cutaneous_survival_OS)], 
  HR.PFS.melanoma = melanoma_survival_PFS[, -ncol(melanoma_survival_PFS)], 
  HR.PFS.cutaneous = cutaneous_survival_PFS[, -ncol(cutaneous_survival_PFS)]
)

Table.S12 <- melanoma_descriptive
  


# Results -----------------------------------------------------------------


if(!dir.exists("../results")) dir.create("../results")
if(!dir.exists("../results/Tables")) dir.create("../results/Tables")
if(!dir.exists("../results/Figures")) dir.create("../results/Figures")

write.table(Table.S1, "../results/Tables/Table-S1.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(Table.S2, "../results/Tables/Table-S2.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(Table.S12, "../results/Tables/Table-S12.txt", col.names = T, row.names = F, sep = "\t", quote = F)
