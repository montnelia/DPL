require("cvAUC")
library("cvAUC")
library("glmnet")
library("pROC")



setwd("C:/Users/ri89por/Desktop/Forschung/SummerRetreat/RCode")
getwd()


getwd()

dataset <- list.files("inst_scRNA")
dataset
pos_dataset <- 2
dataset_desc <- readRDS(paste0("inst_scRNA/", dataset[pos_dataset]))
x <-   dataset_desc[,1:(dim(dataset_desc)[2]-1)]
y <- dataset_desc[,dim(dataset_desc)[2]]
dim(x)
length(y)



#ada_lasso - function(dataset_desc, cluster_function, n_rep, option_tune_gamma, strat){
  n_rep <- 1 # Ziel: 100  (entspricht untigem k) 
  

  list_coef <- list()
  output_acc <- data.frame(matrix(NA, nrow =n_rep+1, ncol = 10))
  colnames(output_acc) <- c("Data", "Lasso", "Ada_Lasso", "TaIl_lasso_DB", "TaIl_lasso_FN", "TaIl_lasso_ANOVA", "TaIl_lasso_sil", "TaIl_lasso_aic", "Ridge", "Enet")
  output_acc_auc <- data.frame(matrix(NA, nrow =n_rep+1, ncol = 10))
  colnames(output_acc_auc) <- c("Data", "Lasso", "Ada_Lasso", "TaIl_lasso_DB", "TaIl_lasso_FN", "TaIl_lasso_ANOVA", "TaIl_lasso_sil", "TaIl_lasso_aic", "Ridge", "Enet")
  
  output_numb_coef <- data.frame(matrix(NA, nrow =n_rep+1, ncol = 10))
  colnames(output_numb_coef) <- c("Data", "Lasso", "Ada_Lasso", "TaIl_lasso_DB", "TaIl_lasso_FN", "TaIl_lasso_ANOVA", "TaIl_lasso_sil", "TaIl_lasso_aic", "Ridge", "Enet")
  
  names(y) <- 1:length(y)
  
  
  set.seed(123)
  for(k in 1:n_rep){
    
    #seeds = sample(1:10000, rep, replace=FALSE)   # Das sind die verschiedenen Versionen, die ich getestet habe und so erschreckend unterschiedliche Performance rauskam...
    #seeds <- 123455:(123455+rep)
    #seed = seeds[k]
    #set.seed(seed)
    
      
      ids             = sample(1:nrow(x)) # not stratified 
      fold            = rep(1:10, length.out = nrow(x))# not stratified 
      
  
    normal_lasso    = c()
    normal_ridge = c()
    normal_enet = c()
    normal_ada_lasso2 = c()
    tail_lasso_db  = c()
    tail_lasso_fn  = c()
    tail_lasso_anova  = c()
    tail_lasso_sil  = c()
    tail_lasso_aic  = c()
    
    normal_lasso_auc    = c()
    normal_enet_auc = c()
    normal_ridge_auc = c()
    normal_ada_lasso2_auc = c()
    tail_lasso_db_auc  = c()
    tail_lasso_fn_auc  = c()
    tail_lasso_anova_auc  = c()
    tail_lasso_sil_auc  = c()
    tail_lasso_aic_auc  = c()
    
    non_zero_normal = c()
    non_zero_ridge = c()
    non_zero_enet = c()
    non_zero_ada2   = c()
    non_zero_tail_lasso_db = c()
    non_zero_tail_lasso_fn = c()
    non_zero_tail_lasso_anova = c()
    non_zero_tail_lasso_sil = c()
    non_zero_tail_lasso_aic  = c()
    
    
   

    

    
    
    for(i in 1:1){ #i in 1:10
  
      rescaled_lasso_nested  = c()
 
        ids_train = ids[fold != i ] 
  

        

      x_train <- x[ids_train,]
      x_train <- as.matrix(x_train)
      
      x_test <- x[-ids_train,]
      x_test <- as.matrix(x_test)
      
      y_train <- as.factor(y[ids_train])
      y_test <- as.factor(y[-ids_train])
      
      ####################
      ### Normal Lasso ###
      ####################
    
      cv_res      = cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)
      pred_lasso  = predict(cv_res, x_test, s = "lambda.min", type="class", alpha =1)
      normal_lasso[i]    = 1 - mean(pred_lasso == y_test)
      non_zero_normal[i]   = sum(cv_res$glmnet.fit$beta[,cv_res$lambda == cv_res$lambda.min] != 0)
      pred_lasso  = predict(cv_res, x_test, s = "lambda.min", type="response", alpha =1)
      normal_lasso_auc[i] <- pROC::auc(y_test, pred_lasso)
      
      
      ### ### Elastic Net ### 
      cv_res      = cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5)
      pred_lasso  = predict(cv_res, x_test, s = "lambda.min", type="class", alpha =0.5)
      normal_enet[i]    = 1 - mean(pred_lasso == y_test)
      non_zero_enet[i]   = sum(cv_res$glmnet.fit$beta[,cv_res$lambda == cv_res$lambda.min] != 0)
      pred_lasso  = predict(cv_res, x_test, s = "lambda.min", type="response", alpha =0.5)
      normal_enet_auc[i] <- pROC::aucauc(y_test, pred_lasso)
      
      ### Ridge ### 
      cv_res      = cv.glmnet(x_train, y_train, family = "binomial", alpha = 0)
      pred_ridge  = predict(cv_res, x_test, s = "lambda.min", type="class", alpha =0)
      normal_ridge[i]    = 1 - mean(pred_lasso == y_test)
      non_zero_ridge[i]   = sum(cv_res$glmnet.fit$beta[,cv_res$lambda == cv_res$lambda.min] != 0)
      pred_lasso  = predict(cv_res, x_test, s = "lambda.min", type="response", alpha =0)
      normal_ridge_auc[i] <- pROC::aucauc(y_test, pred_lasso)
      
      ######################
      ### Adaptive Lasso ###
      ######################
      ridge1_cv <- cv.glmnet(x = x_train, y = y_train,family = "binomial", alpha = 0)
      best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
      
      
      alasso1_cv <- cv.glmnet(x = x_train, y = y_train, family = "binomial", nfold = 10,
                              alpha = 1, penalty.factor = 1 / abs(best_ridge_coef),
                              ## prevalidated array is returned
                              keep = TRUE)
      
      pred_adalasso  = predict(alasso1_cv, x_test, s = "lambda.min", type="class", alpha = 1)
      normal_ada_lasso2[i]    = 1 - mean(pred_adalasso == y_test)
      non_zero_ada2[i] <- sum(alasso1_cv$glmnet.fit$beta[, alasso1_cv$lambda == alasso1_cv$lambda.min] !=0)
      
      pred_lasso  = predict(alasso1_cv, x_test, s = "lambda.min", type="response", alpha =1)
      normal_ada_lasso2_auc[i] <- pROC::aucauc(y_test, pred_lasso)
      
      
    
      ##################
      ### TaIl Lasso ###
      ##################
      ### DB ### 
      cluster_index        = get_DB_index(x[ids_train,], as.integer(y[ids_train]))
      cv_res_weighted      = cv.glmnet(x_train, y_train, family = "binomial", penalty.factor = cluster_index, alpha = 1)
      pred_lasso_weighted  = predict(cv_res_weighted, x_test, s = "lambda.min", type="class", alpha = 1)
      tail_lasso_db[i]    = 1 - mean(pred_lasso_weighted == y[-ids_train])
      non_zero_tail_lasso_db[i]          = sum(cv_res_weighted$glmnet.fit$beta[, cv_res_weighted$lambda == cv_res_weighted$lambda.min] !=0)
      pred_lasso  = predict(cv_res_weighted, x_test, s = "lambda.min", type="response", alpha =1)
      tail_lasso_db_auc[i] <- pROC::aucauc(y_test, pred_lasso)
      
      
      
      ### FN ### 
      cluster_index        = log(1+get_ANOVA_index(x_train, (y_train)))
      
      #cv_res_weighted      = cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, penalty.factor = cluster_index, lambda = seq(0.001,0.005, by =0.0001))
      cv_res_weighted = NA
      while(is.na(cv_res_weighted)){
        tryCatch( { cv_res_weighted <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, penalty.factor = cluster_index)}
                  , error = function(e) {cv_res_weighted <- NA})
      }
      pred_lasso_weighted  = predict(cv_res_weighted, x_test, s = "lambda.min", type="class", alpha = 1)
      tail_lasso_fn[i]    = 1 - mean(pred_lasso_weighted == y[-ids_train])
      non_zero_tail_lasso_fn[i]          = sum(cv_res_weighted$glmnet.fit$beta[, cv_res_weighted$lambda == cv_res_weighted$lambda.min] !=0)
      
      pred_lasso  = predict(cv_res_weighted, x_test, s = "lambda.min", type="response", alpha =1)
      tail_lasso_fn_auc[i] <- pROC::aucauc(y_test, pred_lasso)
      
      
      
      
      ### ANOVA ###
      cluster_index        = log(1+get_ANOVA_index2(x_train, (y_train)))
      
      #cv_res_weighted      = cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, penalty.factor = cluster_index, lambda = seq(0.001,0.005, by =0.0001))
      cv_res_weighted = NA
      while(is.na(cv_res_weighted)){
        tryCatch( { cv_res_weighted <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, penalty.factor = cluster_index)}
                  , error = function(e) {cv_res_weighted <- NA})
      }
      
      pred_lasso_weighted  = predict(cv_res_weighted, x_test, s = "lambda.min", type="class", alpha = 1)
      tail_lasso_anova[i]    = 1 - mean(pred_lasso_weighted == y[-ids_train])
      non_zero_tail_lasso_anova[i]          = sum(cv_res_weighted$glmnet.fit$beta[, cv_res_weighted$lambda == cv_res_weighted$lambda.min] !=0)
      
      pred_lasso  = predict(cv_res_weighted, x_test, s = "lambda.min", type="response", alpha =1)
      tail_lasso_anova_auc[i] <- pROC::aucauc(y_test, pred_lasso)
      
      
      ### Sil ###
      cluster_index        = get_silhouette_index(x[ids_train,], as.integer(y[ids_train]))
      
      cv_res_weighted      = cv.glmnet(x_train, y_train, family = "binomial", penalty.factor = cluster_index, alpha = 1)
      pred_lasso_weighted  = predict(cv_res_weighted, x_test, s = "lambda.min", type="class", alpha = 1)
      tail_lasso_sil[i]    = 1 - mean(pred_lasso_weighted == y[-ids_train])
      non_zero_tail_lasso_sil[i]          = sum(cv_res_weighted$glmnet.fit$beta[, cv_res_weighted$lambda == cv_res_weighted$lambda.min] !=0)
      
      pred_lasso  = predict(cv_res_weighted, x_test, s = "lambda.min", type="response", alpha =1)
      tail_lasso_sil_auc[i] <- pROC::aucauc(y_test, pred_lasso)   
      
      ### AIC ###
      cluster_index        = get_aic_index_binom(x[ids_train,], as.integer(y[ids_train]))
      
      cv_res_weighted      = cv.glmnet(x_train, y_train, family = "binomial", penalty.factor = cluster_index, alpha = 1)
      pred_lasso_weighted  = predict(cv_res_weighted, x_test, s = "lambda.min", type="class", alpha = 1)
      tail_lasso_aic[i]    = 1 - mean(pred_lasso_weighted == y[-ids_train])
      non_zero_tail_lasso_aic[i]          = sum(cv_res_weighted$glmnet.fit$beta[, cv_res_weighted$lambda == cv_res_weighted$lambda.min] !=0)
      
      pred_lasso  = predict(cv_res_weighted, x_test, s = "lambda.min", type="response", alpha =1)
      tail_lasso_aic_auc[i] <- pROC::aucauc(y_test, pred_lasso)   
      
      
      
      
      
      print(paste0(" i =", i))
      
      
    }
    
      

    
    output_acc[k,1] = dataset[pos_dataset]
    output_acc[k,2] = mean(normal_lasso, na.rm = TRUE)
    output_acc[k,3] = mean(normal_ada_lasso2, na.rm = TRUE)
    output_acc[k,4] = mean(tail_lasso_db, na.rm = TRUE)
    output_acc[k,5] = mean(tail_lasso_fn, na.rm = TRUE)
    output_acc[k,6] = mean(tail_lasso_anova, na.rm = TRUE)
    output_acc[k,7] = mean(tail_lasso_sil, na.rm = TRUE)
    output_acc[k,8] = mean(tail_lasso_aic, na.rm = TRUE)
    output_acc[k,9] = mean(normal_ridge, na.rm = TRUE)
    output_acc[k,10] = mean(normal_enet, na.rm = TRUE)
    
    output_acc_auc[k,1] = dataset[pos_dataset]
    output_acc_auc[k,2] = mean(normal_lasso_auc, na.rm = TRUE)
    output_acc_auc[k,3] = mean(normal_ada_lasso2_auc, na.rm = TRUE)
    output_acc_auc[k,4] = mean(tail_lasso_db_auc, na.rm = TRUE)
    output_acc_auc[k,5] = mean(tail_lasso_fn_auc, na.rm = TRUE)
    output_acc_auc[k,6] = mean(tail_lasso_anova_auc, na.rm = TRUE)
    output_acc_auc[k,7] = mean(tail_lasso_sil_auc, na.rm = TRUE)
    output_acc_auc[k,8] = mean(tail_lasso_aic_auc, na.rm = TRUE)
    output_acc_auc[k,9] = mean(normal_ridge_auc, na.rm = TRUE)
    output_acc_auc[k,10] = mean(normal_enet_auc, na.rm = TRUE)
    
    
    output_numb_coef[k,1] = dataset[pos_dataset]
    output_numb_coef[k,2] = mean(non_zero_normal, na.rm = TRUE)
    output_numb_coef[k,3] = mean(non_zero_ada2, na.rm = TRUE)
    output_numb_coef[k,4] = mean(non_zero_tail_lasso_db, na.rm = TRUE)
    output_numb_coef[k,5] = mean(non_zero_tail_lasso_fn, na.rm = TRUE)
    output_numb_coef[k,6] = mean(non_zero_tail_lasso_anova, na.rm = TRUE)
    output_numb_coef[k,7] = mean(non_zero_tail_lasso_sil, na.rm = TRUE)
    output_numb_coef[k,8] = mean(non_zero_tail_lasso_aic, na.rm = TRUE)
    output_numb_coef[k,9] = mean(non_zero_ridge, na.rm = TRUE)
    output_numb_coef[k,10] = mean(non_zero_enet, na.rm = TRUE)
    
    
    

    
  }
  output_numb_coef
  output_acc
  output_acc_auc
  
  saveRDS(output_numb_coef, paste0("Results/Numb_coeff/", dataset[pos_dataset]))
  saveRDS(output_acc_auc, paste0("Results/AUC/", dataset[pos_dataset]))
  saveRDS(output_acc, paste0("Results/Misclassification/", dataset[pos_dataset]))
#}
  
  
  
  