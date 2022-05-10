#########################
### Data Set GSE45719 ###
#########################
#https://www.nature.com/articles/nmeth.4612.pdf NATURE PAPER 
#https://arxiv.org/pdf/1802.09381.pdf Binomial
#https://hal.archives-ouvertes.fr/hal-01716704v2/document Multinomial


#################
### Read data ###
#################
setwd("C:/Users/ri89por/Desktop/Github/TargetAdaptiveLASSO/RCode")
GSE45719 <- readRDS("dataset/GSE45719.rds")
(GSE45719_gene <- experiments(GSE45719)[["gene"]])
data <- t(assays(GSE45719_gene)[["TPM"]])
dim(data)
str(GSE45719)
y_org <- (GSE45719@colData$description)
y_org
length(y_org)
y <- as.factor(c(rep("16-cell stage embryo", 50), rep("4-cell stage embryo",14), 
                 rep("8-cell stage embryo",28), rep("Liver cell", 5),
                 rep("2-cell stage embryo", 8), rep("Early 2-cell stage embryo", 8),
                 rep("Early blastocyst", 43), rep("Late 2-cell stage embryo",10), 
                 rep("Late blastocyst", 30), rep("Mid 2-cell stage embryo", 12),
                 rep("Mid blastocyst", 60), rep("Zygote", 4), rep("8-cell stage embryo", 9),
                 rep("fibroblast", 10)))
#gsub( " . .*$", "", y_org )
y_prep <- y
y_prep
levels(y) <- 1:(length(levels(y)))
x <- apply(data[,1:(dim(data)[2])], 2, as.numeric)
data_aim <- data.frame(x,y)
saveRDS(data_aim, "inst_scRNA/GSE45719.rds")


################
### binomial ###
################
### mid blastocyst and 16 cell ###
binom <- c(which(y_prep=="Mid blastocyst" | y_prep=="16-cell stage embryo"))
length(binom)
x_mid <- x[binom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #10851 for mid blastocyst and 16 cell stage embryo
#freq = apply(x_mid, 2, function(x)mean(round(x) > 1))
#mean(freq > 0.25)

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #13045
y <- as.factor(c(y_prep[which(y_prep=="Mid blastocyst" | y_prep=="16-cell stage embryo")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="Mid blastocyst" | y_prep=="16-cell stage embryo"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/GSE45719_mid_16.rds")



### 8 cell and 16 cell ###
binom <- c(which(y_prep=="8-cell stage embryo" | y_prep=="16-cell stage embryo"))
length(binom)
x_mid <- x[binom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #10626 for 8 cell stage embryo and 16 cell stage embryo
#freq = apply(x_mid, 2, function(x)mean(round(x) > 1))
#mean(freq > 0.25)

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #10626
y <- as.factor(c(y_prep[which(y_prep=="8-cell stage embryo" | y_prep=="16-cell stage embryo")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="8-cell stage embryo" | y_prep=="16-cell stage embryo"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim) # 87 10627
saveRDS(data_aim, "inst_scRNA/GSE45719_8_16.rds")




###################
### multinomial ###
###################

multinom <- c(which(y_prep=="Mid blastocyst" |y_prep=="8-cell stage embryo" | y_prep=="16-cell stage embryo"))
length(multinom)
x_mid <- x[multinom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #11065 

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #11065
y <- as.factor(c(y_prep[which(y_prep=="Mid blastocyst" |y_prep=="8-cell stage embryo" | y_prep=="16-cell stage embryo")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="Mid blastocyst" |y_prep=="8-cell stage embryo" | y_prep=="16-cell stage embryo"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/GSE45719_multi.rds")
