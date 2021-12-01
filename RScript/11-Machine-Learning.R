###################################################
## 
## Machine Learning analysis
##
## - IMK-140 genes with lasso 
## - IMK-RF random forest for gene selection and lasso for training
##
## Manuscript --> Figure 6, Table S11
##
###################################################


library(caret)
library(ranger)
library(glmnet)
library(parallel)

lasso_grid <- expand.grid(
  alpha = 1,
  lambda = seq(0.0001, 0.1, length = 20)
)


rf_grid <- expand.grid(mtry = c(4, 5),
                       splitrule = c("gini"),
                       min.node.size = c(3, 5, 10))

myControl <- trainControl(
  method = "repeatedcv", number = 7, repeats=10,
  summaryFunction = twoClassSummary,
  classProbs = TRUE # Super important!
)



lasso_grid <- expand.grid(
  alpha = 1,
  lambda = seq(0.0001, 0.1, length = 20)
)

myControl <- trainControl(
  method = "repeatedcv", number = 10, repeats=20,
  summaryFunction = twoClassSummary,
  classProbs = TRUE # Super important!
)

# Number of replicates
n <- 500

norm.counts <- read.delim("../data/norm_counts_cutaneous.txt")

# IMK-140
datos <- norm.counts
rownames(datos) <- 1:nrow(datos)
system.time(cv_imk <- mclapply(1:n, cross.validation, mc.cores = detectCores()-1))

# IMK-RF
datos <- norm.counts
rownames(datos) <- 1:nrow(datos)
system.time(cv_imk_rf <- mclapply(1:n, cross.validation_rf, mc.cores = detectCores()-1))

