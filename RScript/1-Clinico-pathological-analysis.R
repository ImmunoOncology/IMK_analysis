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





