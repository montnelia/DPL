
library(clusterSim)
library(ggplot2)
library(cluster)
library(glmnet)
library(clValid)
library(DescTools)

library(gridExtra)
library(lattice)
library("DescTools")



get_DB_index = function(x, y){
  cluster_index = c()
  x_sc = scale(x)
  for(k in 1:ncol(x)){
    cluster_index[k] = index.DB(x_sc[,k], y)$DB
  }
  cluster_index
  }

#db <-  get_DB_index(x,y)


get_silhouette_index = function(x, y){
  cluster_index = c()
  x_sc = scale(x)
  for(k in 1:ncol(x)){
    cluster_index[k] = summary(silhouette(y, dist(x_sc[,k])))$si.summary[4] 
  }
  sil_ind <- abs(1/cluster_index) 
  sil_ind
}
#sil <- get_silhouette_index(x, y)










get_ANOVA_index2 = function(x, y){
  
  
  cluster_index <- c()
  x_sc = scale(x)
  for(k in 1:ncol(x)){ 
    
    data <-data.frame(gene= x_sc[,k], y =y)
    colnames(data) <- c("gene", "y")
    
    anov<- summary(aov(gene ~ y, data = data))
    cluster_index[k] <- anov[[1]]$`F value`[1]
  }
  1/cluster_index
}


#ANOVA <- get_ANOVA_index2(x, y)

#cluster_index        = log(1+get_ANOVA_index2(x_train, (y_train)))

