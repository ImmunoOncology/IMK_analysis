###

do_bivariate_clinical <- function(datos){
  
  tests <- data.frame(Variable = NA, p.value = NA, test = NA)
  
  dummy <- datos[!is.na(datos$Stage), c("Response", "Stage")]
  table(dummy$Response, dummy$Stage)
  p.value <- round(fisher.test(dummy$Response, dummy$Stage)$p.value, 4)
  
  tests <- rbind(tests, data.frame(Variable = "Stage", p.value = p.value, test = "Fisher"))
  tests <- na.omit(tests)
  
  datos$Age
  dummy <- datos[!is.na(datos$Age), c("Response", "Age")]
  shapiro.test(dummy$Age)
  p.value <- round(t.test(dummy$Age~dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "Age", p.value = p.value, test = "T.test"))
  
  datos$Sex
  dummy <- datos[!is.na(datos$Sex), c("Response", "Sex")]
  table(dummy$Response, dummy$Sex)
  p.value <- round(chisq.test(dummy$Sex,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "Sex", p.value = p.value, test = "Chisq"))
  
  
  datos$Hypertension
  dummy <- datos[!is.na(datos$Hypertension), c("Response", "Hypertension")]
  table(dummy$Response, dummy$Hypertension)
  p.value <- round(fisher.test(dummy$Hypertension,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "Hypertension", p.value = p.value, test = "Fisher"))
  
  
  datos$Diabetes
  dummy <- datos[!is.na(datos$Diabetes), c("Response", "Diabetes")]
  table(dummy$Response, dummy$Diabetes)
  p.value <- round(fisher.test(dummy$Diabetes,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "Diabetes", p.value = p.value, test = "Fisher"))
  
  
  datos$DL
  dummy <- datos[!is.na(datos$DL), c("Response", "DL")]
  table(dummy$Response, dummy$DL)
  p.value <- round(fisher.test(dummy$DL,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "DL", p.value = p.value, test = "Fisher"))
  
  
  datos$BRAF
  dummy <- datos[!is.na(datos$BRAF), c("Response", "BRAF")]
  table(dummy$Response, dummy$BRAF)
  p.value <- round(fisher.test(dummy$BRAF,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "BRAF", p.value = p.value, test = "Fisher"))
  
  
  datos$M_Lung
  dummy <- datos[!is.na(datos$M_Lung), c("Response", "M_Lung")]
  table(dummy$Response, dummy$M_Lung)
  p.value <- round(fisher.test(dummy$M_Lung,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "M_Lung", p.value = p.value, test = "Fisher"))
  
  
  datos$M_node
  dummy <- datos[!is.na(datos$M_node), c("Response", "M_node")]
  table(dummy$Response, dummy$M_node)
  p.value <- round(fisher.test(dummy$M_node,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "M_node", p.value = p.value, test = "Fisher"))
  
  
  datos$num_M1
  dummy <- datos[!is.na(datos$num_M1), c("Response", "num_M1")]
  table(dummy$Response, dummy$num_M1)
  p.value <- round(fisher.test(dummy$num_M1,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "num_M1", p.value = p.value, test = "Fisher"))
  
  
  datos$LDH_m1
  dummy <- datos[!is.na(datos$LDH_m1), c("Response", "LDH_m1")]
  shapiro.test(dummy$LDH_m1)
  p.value <- round(wilcox.test(dummy$LDH_m1~dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "LDH_m1", p.value = p.value, test = "wilcox.test"))
  
  
  datos$CNS_metastasis
  dummy <- datos[!is.na(datos$CNS_metastasis), c("Response", "CNS_metastasis")]
  table(dummy$Response, dummy$CNS_metastasis)
  p.value <- round(fisher.test(dummy$CNS_metastasis,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "CNS_metastasis", p.value = p.value, test = "Fisher"))
  
  
  datos$No_Previous_lines
  dummy <- datos[!is.na(datos$No_Previous_lines), c("Response", "No_Previous_lines")]
  table(dummy$Response, dummy$No_Previous_lines)
  p.value <- round(fisher.test(dummy$No_Previous_lines,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "No_Previous_lines", p.value = p.value, test = "Fisher"))
  
  
  datos$LDH_treatment_previous_to_IT
  dummy <- datos[!is.na(datos$LDH_treatment_previous_to_IT), c("Response", "LDH_treatment_previous_to_IT")]
  shapiro.test(dummy$LDH_treatment_previous_to_IT)
  p.value <- round(wilcox.test(dummy$LDH_treatment_previous_to_IT~dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "LDH_treatment_previous_to_IT", p.value = p.value, test = "wilcox.test"))
  
  
  datos$Neutrophiles_treatment_previous_to_IT
  dummy <- datos[!is.na(datos$Neutrophiles_treatment_previous_to_IT), c("Response", "Neutrophiles_treatment_previous_to_IT")]
  shapiro.test(dummy$Neutrophiles_treatment_previous_to_IT)
  p.value <- round(wilcox.test(dummy$Neutrophiles_treatment_previous_to_IT~dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "Neutrophiles_treatment_previous_to_IT", p.value = p.value, test = "wilcox.test"))
  
  
  datos$Lymphocytes_treatment_previous_to_IT
  dummy <- datos[!is.na(datos$Lymphocytes_treatment_previous_to_IT), c("Response", "Lymphocytes_treatment_previous_to_IT")]
  shapiro.test(dummy$Lymphocytes_treatment_previous_to_IT)
  p.value <- round(t.test(dummy$Lymphocytes_treatment_previous_to_IT~dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "Lymphocytes_treatment_previous_to_IT", p.value = p.value, test = "T.test"))
  
  
  
  datos$Platelets_treatment_previous_to_IT
  dummy <- datos[!is.na(datos$Platelets_treatment_previous_to_IT), c("Response", "Platelets_treatment_previous_to_IT")]
  shapiro.test(dummy$Platelets_treatment_previous_to_IT)
  p.value <- round(t.test(dummy$Platelets_treatment_previous_to_IT~dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "Platelets_treatment_previous_to_IT", p.value = p.value, test = "T.test"))
  
  
  datos$Toxicity_IT
  dummy <- datos[!is.na(datos$Toxicity_IT), c("Response", "Toxicity_IT")]
  table(dummy$Response, dummy$Toxicity_IT)
  p.value <- round(fisher.test(dummy$Toxicity_IT,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "Toxicity_IT", p.value = p.value, test = "Fisher"))
  
  datos$Maximum_toxicity_grade
  dummy <- datos[!is.na(datos$Maximum_toxicity_grade), c("Response", "Maximum_toxicity_grade")]
  table(dummy$Response, dummy$Maximum_toxicity_grade)
  p.value <- round(fisher.test(dummy$Maximum_toxicity_grade,dummy$Response)$p.value, 4)
  tests <- rbind(tests, data.frame(Variable = "Maximum_toxicity_grade", p.value = p.value, test = "Fisher"))
  
  tests <- tests[order(tests$p.value, decreasing = F), ]
  return(tests)
  
}


do_descriptive <- function(datos, vars){

  result <- data.frame(my.vars = vars)
  rownames(result) <- result$my.vars
  
  for(var in vars){
    if(is.numeric(datos[, var])){
      result[var, "Value"] <- round(median(datos[, var], na.rm = T), 2)
      result[var, "Good.value"] <- round(median(datos[datos$Response=="Good", var], na.rm = T), 2)
      result[var, "Bad.value"] <- round(median(datos[datos$Response=="Bad", var], na.rm = T), 2)
      result[var, "Range"] <- paste(round(range(datos[, var], na.rm = T), 2), collapse = " - ")
      result[var, "Good.range"] <- paste(round(range(datos[datos$Response=="Good", var], na.rm = T), 2), collapse = " - ")
      result[var, "Bad.range"] <- paste(round(range(datos[datos$Response=="Bad", var], na.rm = T), 2), collapse = " - ")
    }else{
      tab <- table(datos[, var])
      for(i in rownames(tab)){
        result <- rbind(result, NA)
        result$my.vars[nrow(result)] <- paste(var, i, sep = "_")
        rownames(result)[nrow(result)] <- paste(var, i, sep = "_")
        result[paste(var, i, sep = "_"), "Value"] <- tab[i]
        result[paste(var, i, sep = "_"), "Range"] <- round(tab[i]/sum(tab)*100, 2)
        result[paste(var, i, sep = "_"), "sep"]<-T
      }
      tab <- table(datos[datos$Response=="Good", var])
      for(i in rownames(tab)){
        result[paste(var, i, sep = "_"), "Good.value"] <- tab[i]
        result[paste(var, i, sep = "_"), "Good.range"] <- round(tab[i]/sum(tab)*100, 2)
      }
      tab <- table(datos[datos$Response=="Bad", var])
      for(i in rownames(tab)){
        result[paste(var, i, sep = "_"), "Bad.value"] <- tab[i]
        result[paste(var, i, sep = "_"), "Bad.range"] <- round(tab[i]/sum(tab)*100, 2)
      }
      
    }
    
    pct.na <- sum(is.na(datos[, var]))
    result[var, "NA pct"] <- paste(pct.na, " (", round(pct.na/length(datos[, var]), 2), "%)", sep = "")
    pct.na <- sum(is.na(datos[datos$Response=="Good", var]))
    result[var, "Good NA pct"] <- paste(pct.na, " (", round(pct.na/length(datos[, var]), 2), "%)", sep = "")
    pct.na <- sum(is.na(datos[datos$Response=="Bad", var]))
    result[var, "Bad NA pct"] <- paste(pct.na, " (", round(pct.na/length(datos[, var]), 2), "%)", sep = "")
    
  }
  
  result <- result[order(result$my.vars), ]
  
  idt <- grep("-", result$Range)
  result$Value[idt] <- paste(result$Value[idt], " [",
                             result$Range[idt], "]", sep = "")
  result$Value[-idt] <- paste(result$Value[-idt], " (",
                              result$Range[-idt], " %)", sep = "")
  result$Good.value[idt] <- paste(result$Good.value[idt], " [",
                                  result$Good.range[idt], "]", sep = "")
  result$Good.value[-idt] <- paste(result$Good.value[-idt], " (",
                                   result$Good.range[-idt], " %)", sep = "")
  result$Bad.value[idt] <- paste(result$Bad.value[idt], " [",
                                 result$Bad.range[idt], "]", sep = "")
  result$Bad.value[-idt] <- paste(result$Bad.value[-idt], " (",
                                  result$Bad.range[-idt], " %)", sep = "")
  
  result$Bad.range <- NULL
  result$Good.range <- NULL
  result$Range <- NULL
  
  return(result)
}

do_survival <- function(datos, vars, time, event){
  
  univ_formulas <- sapply(vars, function(x) as.formula(paste('Surv(', time, ',', event, ')~', x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = datos)})
  uni_survival <- forestmodel::forest_model(model_list = univ_models,covariates = vars,merge_models = T, return_data=T)
  uni_survival <- uni_survival$plot_data$forest_data
  uni_survival$variable <- gsub("[._]", " ", uni_survival$variable)
  uni_survival$variable[is.na(uni_survival$variable)] <- ""
  uni_survival$coef <- "Reference"
  idt <- !uni_survival$reference
  uni_survival$estimate <- round(exp(uni_survival$estimate), 4)
  uni_survival$conf.low <- round(exp(uni_survival$conf.low), 4)
  uni_survival$conf.high <- round(exp(uni_survival$conf.high), 4)
  uni_survival$coef[idt] <- paste0(uni_survival$estimate[idt], " ", "(", uni_survival$conf.low[idt], " - ", uni_survival$conf.high[idt], ")")
  uni_survival <- uni_survival[, c("variable", "level", "n", "coef", "p.value")]
  uni_survival$p.value <- round(uni_survival$p.value, 4)
  uni_logrank <- lapply(univ_formulas, function(x){surv_pvalue(surv_fit(x, data = datos))$pval})
  uni_survival$Logrank <- NA
  uni_survival$Logrank[match(names(uni_logrank), names(uni_survival$variable))] <- unlist(uni_logrank)
  return(uni_survival)
}

do_ggsurvival_gene <- function(df_surv, time, event, my.col  = c("Low"="firebrick2", "Med"="mediumseagreen", "High"="skyblue"), xlab = "OS_event in days"
                                                               , ylab="Overall survival probability", output_dir=NULL){
  
  res <- list()
  for(i in colnames(event_rna)){
    
    if(length(table(df_surv[,i]))>1){
      
      a=survdiff(Surv(df_surv[,time], df_surv[,event]) ~ df_surv[,i],data=df_surv)
      res[[i]] <- 1 - pchisq(a$chisq, length(a$n) - 1)
      
      if(!is.null(output_dir)){
        if(res[[i]]<0.05){
          df_surv$var <- factor(df_surv[, i], levels = c("low", "med", "high"))
          levels(df_surv$var) <- c("Low", "Med", "High")
          
          df_surv$my_time <- df_surv[, time]
          df_surv$my_event <- df_surv[, event]
          
          p <- ggsurvplot(
            survfit(Surv(my_time, my_event) ~ var, data = df_surv),                 
            conf.int = F,         
            size=1,#0.7,                    
            pval=TRUE, 
            palette = my.col[levels(df_surv$var)[levels(df_surv$var)%in%unique(df_surv$var)]],
            pval.method=TRUE,   break.x.by = 500,    
            xlab=xlab,
            ylab=ylab,
            ylim=c(0,1),
            surv.scale="percent",
            legend.labs = levels(df_surv$var)[levels(df_surv$var)%in%unique(df_surv$var)],
            tables.col="strata",
            risk.table=T, 
            risk.table.col = "strata",
            risk.table.y.text = T,
            legend = c(.7,.8),
            tables.y.text = T, 
            legend.title="", ggtheme = custom_theme(), title = i)
        ggsave(paste0(output_dir, i, ".pdf"), plot = print(p), width = 6.5, height = 7, onefile=F)
        }
      }
    }
  }
  
  return(res)
}


cross.validation <- function(i){
  
  result <- list()
  result[["prob"]] <- list()
  result[["resp"]] <- list()
  result[["coef"]] <- list()
  result[["lambda"]] <- list()
  
  idt <- createFolds(datos$response, k = 10)
  cnt <- 1
  
  for(j in idt){
    
    train <- datos[-j, ]
    test <- datos[j, ]
    
    model.lasso <- train(response ~ ., train, method = "glmnet", 
                         trControl = myControl, tuneGrid = lasso_grid)
    
    model <- glmnet::glmnet(model.matrix(response~., train), y = train$response, 
                            family = "binomial", alpha = model.lasso$bestTune$alpha, lambda = model.lasso$bestTune$lambda)
    
    result[["prob"]][[cnt]] <- predict(model, newx = model.matrix(response~., test), type = "response")
    result[["resp"]][[cnt]] <- test$response
    result[["coef"]][[cnt]] <- as.matrix(coef(model))
    result[["lambda"]][[cnt]] <- model.lasso$bestTune$lambda
    cnt <- cnt + 1
    message(cnt)
  }
  
  return(result)
}


cross.validation <- function(i){
  
  result <- list()
  result[["prob"]] <- list()
  result[["resp"]] <- list()
  result[["coef"]] <- list()
  result[["lambda"]] <- list()
  
  idt <- createFolds(datos$response, k = 10)
  cnt <- 1
  
  for(j in idt){
    
    train <- datos[-j, ]
    test <- datos[j, ]
    
    model.lasso <- train(response ~ ., train, method = "glmnet", 
                         trControl = myControl, tuneGrid = lasso_grid)
    
    model <- glmnet::glmnet(model.matrix(response~., train), y = train$response, 
                            family = "binomial", alpha = model.lasso$bestTune$alpha, lambda = model.lasso$bestTune$lambda)
    
    result[["prob"]][[cnt]] <- predict(model, newx = model.matrix(response~., test), type = "response")
    result[["resp"]][[cnt]] <- test$response
    result[["coef"]][[cnt]] <- as.matrix(coef(model))
    result[["lambda"]][[cnt]] <- model.lasso$bestTune$lambda
    cnt <- cnt + 1
    message(cnt)
  }
  
  return(result)
}

cross.validation_rf <- function(i){
  
  result <- list()
  result[["prob"]] <- list()
  result[["resp"]] <- list()
  result[["coef"]] <- list()
  result[["lambda"]] <- list()
  
  idt <- createFolds(datos$response, k = 10)
  cnt <- 1
  
  for(j in idt){
    
    train <- datos[-j, ]
    test <- datos[j, ]
    
    rf_fit <- train(response ~ ., 
                    data = train, 
                    method = "ranger",
                    trControl = myControl,
                    importance = "permutation",
                    tuneGrid = rf_grid)
    
    train <- train[, c(names(varImp(rf_fit)$importance)[1:20], "response")]
    
    model.lasso <- train(response ~ ., train, method = "glmnet", 
                         trControl = myControl, tuneGrid = lasso_grid)
    
    model <- glmnet::glmnet(model.matrix(response~., train), y = train$response, 
                            family = "binomial", alpha = model.lasso$bestTune$alpha, lambda = model.lasso$bestTune$lambda)
    
    result[["prob"]][[cnt]] <- predict(model, newx = model.matrix(response~., test), type = "response")
    result[["resp"]][[cnt]] <- test$response
    result[["coef"]][[cnt]] <- as.matrix(coef(model))
    result[["lambda"]][[cnt]] <- model.lasso$bestTune$lambda
    cnt <- cnt + 1
    message(cnt)
  }
  
  return(result)
}


