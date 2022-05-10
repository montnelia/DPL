require("cvAUC")
library("cvAUC")
library("glmnet")
library("pROC")
library(HandTill2001) # Multinomial AUC



setwd("C:/Users/ri89por/Desktop/TargetAdaptiveLASSO/RCode")
getwd()

dataset <- list.files("inst_scRNA")
dataset
pos_dataset <- 3
dataset_desc <- readRDS(paste0("inst_scRNA/", dataset[pos_dataset]))
x <-   dataset_desc[,1:(dim(dataset_desc)[2]-1)]
y <- dataset_desc[,dim(dataset_desc)[2]]



#ada_lasso - function(dataset_desc, cluster_function, n_rep, option_tune_gamma, strat){
  n_rep <- 5# Ziel: 100  (entspricht untigem k) 
  

  list_coef <- list()
  output_acc <- data.frame(matrix(NA, nrow =n_rep+1, ncol = 9))
  colnames(output_acc) <- c("Data", "Lasso", "Ada_Lasso", "TaIl_lasso_DB", "TaIl_lasso_ANOVA", "TaIl_lasso_sil", "TaIl_lasso_aic", "Ridge", "Enet")
  output_acc_auc <- data.frame(matrix(NA, nrow =n_rep+1, ncol = 9))
  colnames(output_acc_auc) <- c("Data", "Lasso", "Ada_Lasso", "TaIl_lasso_DB",  "TaIl_lasso_ANOVA", "TaIl_lasso_sil", "TaIl_lasso_aic", "Ridge", "Enet")
  
  output_numb_coef <- data.frame(matrix(NA, nrow =n_rep+1, ncol = 9))
  colnames(output_numb_coef) <- c("Data", "Lasso", "Ada_Lasso", "TaIl_lasso_DB",  "TaIl_lasso_ANOVA", "TaIl_lasso_sil", "TaIl_lasso_aic", "Ridge", "Enet")
  
  output_comptime <- data.frame(matrix(NA, nrow =n_rep+1, ncol = 9))
  colnames(output_comptime) <- c("Data", "Lasso", "Ada_Lasso", "TaIl_lasso_DB",  "TaIl_lasso_ANOVA", "TaIl_lasso_sil", "TaIl_lasso_aic", "Ridge", "Enet")
  

  names(y) <- 1:length(y)
  
  output_sel_genes_lasso = list()
  output_sel_genes_ada2   = list()
  output_sel_genes_tail_lasso_db = list()
  output_sel_genes_tail_lasso_fn = list()
  output_sel_genes_tail_lasso_anova = list()
  output_sel_genes_tail_lasso_sil = list()
  output_sel_genes_tail_lasso_aic = list()
  output_sel_genes_ridge = list()
  output_sel_genes_enet = list()
  
  
  
  
  names(y) <- 1:length(y)
  
  
  set.seed(123)
  #set.seed(170192)
  for(r in 1:n_rep){
    
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
    
    time_lasso = c()
    time_ridge = c()
    time_enet = c()
    time_ada2   = c()
    time_tail_lasso_db = c()
    time_tail_lasso_fn = c()
    time_tail_lasso_anova = c()
    time_tail_lasso_sil = c()
    time_tail_lasso_aic  = c()
    
    sel_genes_lasso = vector(mode = "list", length = 10)
    sel_genes_ridge = vector(mode = "list", length = 10)
    sel_genes_enet = vector(mode = "list", length = 10)
    sel_genes_ada2   = vector(mode = "list", length = 10)
    sel_genes_tail_lasso_db = vector(mode = "list", length = 10)
    sel_genes_tail_lasso_fn = vector(mode = "list", length = 10)
    sel_genes_tail_lasso_anova =vector(mode = "list", length = 10)
    sel_genes_tail_lasso_sil = vector(mode = "list", length = 10)
    sel_genes_tail_lasso_aic = vector(mode = "list", length = 10)
    

    
    
    for(i in 2:10){ #i in 1:10
  
      rescaled_lasso_nested  = c()
 
        ids_train = ids[fold != i ] 
  

        

      x_train <- x[ids_train,]
      x_train <- as.matrix(x_train)
      
      x_test <- x[-ids_train,]
      x_test <- as.matrix(x_test)
      
      y_train <- as.factor(y[ids_train])
      y_test <- as.factor(y[-ids_train])
      
      
      #######################
      ### Classical LASSO ###
      #######################
      start_lasso <- Sys.time()
      cv_res      = cv.glmnet(x_train, y_train, family = "multinomial", alpha = 1) #included CV
      pred_lasso  = predict(cv_res, x_test, s = "lambda.min", type="class", alpha =1)
      normal_lasso[i]    = 1 - mean(pred_lasso == y_test)
    
      pred_lasso  = predict(cv_res, x_test, s = "lambda.min", type="response", alpha =1)
      normal_lasso_auc[i] <- HandTill2001::auc(multcap(
        response = y_test,
        predicted = as.matrix(pred_lasso[,,1])
      ))
      
      
      coef_names <- c()
      for(l in 1:length(cv_res$glmnet.fit$beta)){
       
        coef_tmp <- cv_res$glmnet.fit$beta[[l]][,which(cv_res$lambda == cv_res$lambda.min)]
        coef_tmp
        coef_names <- c(coef_names, names(which(coef_tmp != 0)))
      }
      
      non_zero_normal[i] <- length(unique(coef_names))
     
      
      sel_genes_lasso[[i]] <- unique(coef_names)
    
      end_lasso <- Sys.time()
      time_lasso[i] <-  end_lasso - start_lasso
      
      
      
      
      ### Elastic Net ### 
      start_lasso <- Sys.time()
      cv_res      = cv.glmnet(x_train, y_train, family = "multinomial", alpha = 0.5)
      pred_lasso  = predict(cv_res, x_test, s = "lambda.min", type="class", alpha =0.5)
      normal_enet[i]    = 1 - mean(pred_lasso == y_test)
      pred_lasso  = predict(cv_res, x_test, s = "lambda.min", type="response", alpha =0.5)
      normal_enet_auc[i] <- HandTill2001::auc(multcap(
        response = y_test,
        predicted = as.matrix(pred_lasso[,,1])
      ))
      
      coef_names <- c()
      for(l in 1:length(cv_res$glmnet.fit$beta)){
        coef_tmp <- cv_res$glmnet.fit$beta[[l]][,which(cv_res$lambda == cv_res$lambda.min)]
        coef_tmp
        coef_names <- c(coef_names, names(which(coef_tmp != 0)))
      }
      non_zero_enet[i] <- length(unique(coef_names))
      
      
      sel_genes_enet[[i]] <- unique(coef_names)
      
      end_lasso <- Sys.time()
      time_enet[i] <-  end_lasso - start_lasso
      
      ### Ridge ### 
      #start_lasso <- Sys.time()
      #cv_res      = cv.glmnet(x_train, y_train, family = "multinomial", alpha = 0)
      #pred_ridge  = predict(cv_res, x_test, s = "lambda.min", type="class", alpha =0)
      #normal_ridge[i]    = 1 - mean(pred_lasso == y_test)
      #pred_lasso  = predict(cv_res, x_test, s = "lambda.min", type="response", alpha =0)
     
      #normal_ridge_auc[i] <-  HandTill2001::auc(multcap(
       # response = y_test,
        #predicted = as.matrix(pred_lasso[,,1])
      #))
      
      #coef_names <- c()
      #for(l in 1:length(cv_res$glmnet.fit$beta)){
      #  coef_tmp <- cv_res$glmnet.fit$beta[[l]][,which(cv_res$lambda == cv_res$lambda.min)]
       # coef_tmp
      #  coef_names <- c(coef_names, names(which(coef_tmp != 0)))
      #}
      #non_zero_ridge[i] <- length(unique(coef_names))
      
      
      #sel_genes_ridge[[i]] <- unique(coef_names)
      
      
      #end_lasso <- Sys.time()
      #time_ridge[i] <-  end_lasso - start_lasso
      
      
      
      
      ######################
      ### Adaptive Lasso ###
      ######################
      start_lasso <- Sys.time()
      ridge <- cv.glmnet(x = x_train, y_train, alpha = 0, family="multinomial",  type.measure = 'class', grouped=FALSE)
      best_ridge_coef2 <- do.call(cbind, coef(ridge, s = ridge$lambda.min))
      best_ridge_coef2
      best_ridge_weights <- 1 / abs(as.matrix(best_ridge_coef2)[-1,]) # without intercept 
      best_ridge_weights
      
      
      
      
      alasso1_cv <- cv.glmnet(x = x_train, y = y_train, family = "multinomial", nfold = 10,
                             alpha = 1, penalty.factor = (best_ridge_weights),
                            ## prevalidated array is returned
                             keep = TRUE)
      
      pred_adalasso  = predict(alasso1_cv, x_test, s = "lambda.min", type="class", alpha = 1) # type.multinomial = "grouped"
      normal_ada_lasso2[i]    = 1 - mean(pred_adalasso == y_test)
      
      
      pred_lasso  = predict(alasso1_cv, x_test, s = "lambda.min", type="response", alpha =1)
      normal_ada_lasso2_auc[i] <- HandTill2001::auc(multcap(
        response = y_test,
        predicted = as.matrix(pred_lasso[,,1])
      ))
      
      coef_names <- c()
      for(l in 1:length(alasso1_cv$glmnet.fit$beta)){
        
        coef_tmp <- alasso1_cv$glmnet.fit$beta[[l]][,alasso1_cv$lambda == alasso1_cv$lambda.min]
        coef_tmp
        coef_names <- c(coef_names, which(coef_tmp != 0))
      }
      
      non_zero_ada2[i] <- length(unique((coef_names)))

      
      sel_genes_ada2[[i]] = unique(names(coef_names))
    
     
      
      end_lasso <- Sys.time()
      time_ada2[i] <-  end_lasso - start_lasso
      
      
    
      ##################
      ### TaIl Lasso ###
      ##################
      ### DB ### 
      start_lasso <- Sys.time()
      cluster_index        = get_DB_index(x[ids_train,], as.integer(y[ids_train]))
      str(cluster_index)
      cv_res_weighted      = cv.glmnet(x_train, y_train, family = "multinomial", penalty.factor = cluster_index, alpha = 1)
      pred_lasso_weighted  = predict(cv_res_weighted, x_test, s = "lambda.min", type="class", alpha = 1)
      
      tail_lasso_db[i]    = 1 - mean(pred_lasso_weighted == y[-ids_train])
      pred_lasso  = predict(cv_res_weighted, x_test, s = "lambda.min", type="response", alpha =1)
      tail_lasso_db_auc[i] <-  HandTill2001::auc(multcap(
                                              response = y_test,
                                              predicted = as.matrix(pred_lasso[,,1])
                                              ))
      
      
      coef_names <- c()
      for(l in 1:length(cv_res_weighted$glmnet.fit$beta)){
        
        coef_tmp <- cv_res_weighted$glmnet.fit$beta[[l]][,cv_res_weighted$lambda == cv_res_weighted$lambda.min]
        coef_tmp
        coef_names <- c(coef_names, which(coef_tmp != 0))
      }
      
      non_zero_tail_lasso_db[i] <- length(unique(names(coef_names)))
    
      sel_genes_tail_lasso_db[[i]] = unique(names(coef_names))
      
      
      
      end_lasso <- Sys.time()
      time_tail_lasso_db[i] <-  end_lasso - start_lasso
      
      
      ### ANOVA ###
     
      start_lasso <- Sys.time()
      cluster_index        = get_ANOVA_index2(x_train, (y_train))
      str(cluster_index)
      cv_res_weighted      = cv.glmnet(x_train, y_train, family = "multinomial", penalty.factor = cluster_index, alpha = 1)
      pred_lasso_weighted  = predict(cv_res_weighted, x_test, s = "lambda.min", type="class", alpha = 1)
      tail_lasso_anova[i]    = 1 - mean(pred_lasso_weighted == y[-ids_train])
      
      pred_lasso  = predict(cv_res_weighted, x_test, s = "lambda.min", type="response", alpha =1)
      tail_lasso_anova_auc[i] <-  HandTill2001::auc(multcap(
                                                  response = y_test,
                                                  predicted = as.matrix(pred_lasso[,,1])
                                                  ))
      
      coef_names <- c()
      for(l in 1:length(cv_res_weighted$glmnet.fit$beta)){
        
        coef_tmp <- cv_res_weighted$glmnet.fit$beta[[l]][,cv_res_weighted$lambda == cv_res_weighted$lambda.min]
        coef_tmp
        coef_names <- c(coef_names, which(coef_tmp != 0))
      }
      
      non_zero_tail_lasso_anova[i] <- length(unique(names(coef_names)))
      
      sel_genes_tail_lasso_anova[[i]] = unique(names(coef_names))
    
      
      end_lasso <- Sys.time()
      time_tail_lasso_anova[i] <-  end_lasso - start_lasso
      
      
      
      ### Sil ###
      start_lasso <- Sys.time()
      cluster_index        = get_silhouette_index(x[ids_train,], as.integer(y[ids_train]))
      
      cv_res_weighted      = cv.glmnet(x_train, y_train, family = "multinomial", penalty.factor = cluster_index, alpha = 1)
      pred_lasso_weighted  = predict(cv_res_weighted, x_test, s = "lambda.min", type="class", alpha = 1)
      tail_lasso_sil[i]    = 1 - mean(pred_lasso_weighted == y[-ids_train])
      
      pred_lasso  = predict(cv_res_weighted, x_test, s = "lambda.min", type="response", alpha =1)
      tail_lasso_sil_auc[i] <-  HandTill2001::auc(multcap(
                                                response = y_test,
                                                predicted = as.matrix(pred_lasso[,,1])
                                                ))
      
      coef_names <- c()
      for(l in 1:length(cv_res_weighted$glmnet.fit$beta)){
        
        coef_tmp <- cv_res_weighted$glmnet.fit$beta[[l]][,cv_res_weighted$lambda == cv_res_weighted$lambda.min]
        coef_tmp
        coef_names <- c(coef_names, which(coef_tmp != 0))
      }
      
      non_zero_tail_lasso_sil[i] <- length(unique(names(coef_names)))
      
      sel_genes_tail_lasso_sil[[i]] = unique(names(coef_names))
      
    
      
      end_lasso <- Sys.time()
      time_tail_lasso_sil[i] <-  end_lasso - start_lasso
      
      ### AIC ###
      start_lasso <- Sys.time()
      cluster_index        = get_aic_index(x[ids_train,], as.integer(y[ids_train]))
      
      cv_res_weighted      = cv.glmnet(x_train, y_train, family = "multinomial", penalty.factor = cluster_index, alpha = 1)
      pred_lasso_weighted  = predict(cv_res_weighted, x_test, s = "lambda.min", type="class", alpha = 1)
      tail_lasso_aic[i]    = 1 - mean(pred_lasso_weighted == y[-ids_train])
      
      pred_lasso  = predict(cv_res_weighted, x_test, s = "lambda.min", type="response", alpha =1)
      tail_lasso_aic_auc[i] <-  HandTill2001::auc(multcap(
        response = y_test,
        predicted = as.matrix(pred_lasso[,,1])
      ))
      
      
      
      coef_names <- c()
      for(l in 1:length(cv_res_weighted$glmnet.fit$beta)){
        
        coef_tmp <- cv_res_weighted$glmnet.fit$beta[[l]][,cv_res_weighted$lambda == cv_res_weighted$lambda.min]
        coef_tmp
        coef_names <- c(coef_names, which(coef_tmp != 0))
      }
      
      non_zero_tail_lasso_aic[i] <- length(unique(names(coef_names)))
      
      sel_genes_tail_lasso_aic[[i]] = unique(names(coef_names))
      
      end_lasso <- Sys.time()
      time_tail_lasso_aic[i] <-  end_lasso - start_lasso
      
      print(paste0(" i =", i))
      
      
    }
    
      

    
    output_acc[r,1] = dataset[pos_dataset]
    output_acc[r,2] = mean(normal_lasso, na.rm = TRUE)
    output_acc[r,3] = mean(normal_ada_lasso2, na.rm = TRUE)
    output_acc[r,4] = mean(tail_lasso_db, na.rm = TRUE)
    output_acc[r,5] = mean(tail_lasso_anova, na.rm = TRUE)
    output_acc[r,6] = mean(tail_lasso_sil, na.rm = TRUE)
    output_acc[r,7] = mean(tail_lasso_aic, na.rm = TRUE)
    output_acc[r,8] = mean(normal_ridge, na.rm = TRUE)
    output_acc[r,9] = mean(normal_enet, na.rm = TRUE)
    
    
    output_acc_auc[r,1] = dataset[pos_dataset]
    output_acc_auc[r,2] = mean(normal_lasso_auc, na.rm = TRUE)
    output_acc_auc[r,3] = mean(normal_ada_lasso2_auc, na.rm = TRUE)
    output_acc_auc[r,4] = mean(tail_lasso_db_auc, na.rm = TRUE)
    output_acc_auc[r,5] = mean(tail_lasso_anova_auc, na.rm = TRUE)
    output_acc_auc[r,6] = mean(tail_lasso_sil_auc, na.rm = TRUE)
    output_acc_auc[r,7] = mean(tail_lasso_aic_auc, na.rm = TRUE)
    output_acc_auc[r,8] = mean(normal_ridge_auc, na.rm = TRUE)
    output_acc_auc[r,9] = mean(normal_enet_auc, na.rm = TRUE)
    
    
    output_numb_coef[r,1] = dataset[pos_dataset]
    output_numb_coef[r,2] = mean(non_zero_normal, na.rm = TRUE)
    output_numb_coef[r,3] = mean(non_zero_ada2, na.rm = TRUE)
    output_numb_coef[r,4] = mean(non_zero_tail_lasso_db, na.rm = TRUE)
    output_numb_coef[r,5] = mean(non_zero_tail_lasso_anova, na.rm = TRUE)
    output_numb_coef[r,6] = mean(non_zero_tail_lasso_sil, na.rm = TRUE)
    output_numb_coef[r,7] = mean(non_zero_tail_lasso_aic, na.rm = TRUE)
    output_numb_coef[r,8] = mean(non_zero_ridge, na.rm = TRUE)
    output_numb_coef[r,9] = mean(non_zero_enet, na.rm = TRUE)


    output_comptime[r,1] = dataset[pos_dataset]
    output_comptime[r,2] = mean(time_lasso, na.rm = TRUE)
    output_comptime[r,3] = mean(time_ada2, na.rm = TRUE)
    output_comptime[r,4] = mean(time_tail_lasso_db, na.rm = TRUE)
    output_comptime[r,5] = mean(time_tail_lasso_anova, na.rm = TRUE)
    output_comptime[r,6] = mean(time_tail_lasso_sil, na.rm = TRUE)
    output_comptime[r,7] = mean(time_tail_lasso_aic, na.rm = TRUE)
    output_comptime[r,8] = mean(time_ridge, na.rm = TRUE)
    output_comptime[r,9] = mean(time_enet, na.rm = TRUE)
    
    output_sel_genes_lasso[[r]] = unlist(sel_genes_lasso)
    output_sel_genes_ada2[[r]]   = unlist(sel_genes_ada2)
    output_sel_genes_tail_lasso_db[[r]] = unlist(sel_genes_tail_lasso_db)
    output_sel_genes_tail_lasso_fn[[r]] = unlist(sel_genes_tail_lasso_fn)
    output_sel_genes_tail_lasso_anova[[r]] = unlist(sel_genes_tail_lasso_anova)
    output_sel_genes_tail_lasso_sil[[r]] = unlist(sel_genes_tail_lasso_sil)
    output_sel_genes_tail_lasso_aic[[r]] = unlist(sel_genes_tail_lasso_aic)
    output_sel_genes_ridge[[r]] = unlist(sel_genes_ridge)
    output_sel_genes_enet[[r]] = unlist(sel_genes_enet)
    
    
    
    
  }
  output_numb_coef
  output_acc
  output_acc_auc
  
  saveRDS(output_numb_coef, paste0("Results/Numb_coeff_multi/", dataset[pos_dataset]))
  saveRDS(output_acc_auc, paste0("Results/AUC_multi/", dataset[pos_dataset]))
  saveRDS(output_acc, paste0("Results/Misclassification_multi/", dataset[pos_dataset]))
  saveRDS(output_comptime, paste0("Results/Comptime_multi/", dataset[pos_dataset]))
  
  saveRDS(output_sel_genes_lasso, paste0("Results/Sel_genes_multi/Sel_genes_lasso/", dataset[pos_dataset]))
  saveRDS(output_sel_genes_ada2, paste0("Results/Sel_genes_multi/Sel_genes_ada2/", dataset[pos_dataset]))
  saveRDS(output_sel_genes_tail_lasso_db, paste0("Results/Sel_genes_multi/Sel_genes_tail_lasso_db/", dataset[pos_dataset]))
  saveRDS(output_sel_genes_tail_lasso_fn, paste0("Results/Sel_genes_multi/Sel_genes_tail_lasso_fn/", dataset[pos_dataset]))
  saveRDS(output_sel_genes_tail_lasso_anova, paste0("Results/Sel_genes_multi/Sel_genes_tail_lasso_anova/", dataset[pos_dataset]))
  saveRDS(output_sel_genes_tail_lasso_sil, paste0("Results/Sel_genes_multi/Sel_genes_tail_lasso_sil/", dataset[pos_dataset]))
  saveRDS(output_sel_genes_tail_lasso_aic, paste0("Results/Sel_genes_multi/Sel_genes_tail_lasso_aic/", dataset[pos_dataset]))
  #saveRDS(output_sel_genes_ridge, paste0("Results/Sel_genes_multi/Sel_genes_ridge", dataset[pos_dataset]))
  saveRDS(output_sel_genes_enet, paste0("Results/Sel_genes_multi/Sel_genes_tail_lasso_enet", dataset[pos_dataset]))
  
  
  
 
#}
  
  