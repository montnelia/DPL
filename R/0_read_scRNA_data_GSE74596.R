

#########################
### Data Set GSE74596 ###
#########################
GSE74596 <- readRDS("dataset/GSE74596.rds")
GSE74596_tx <- experiments(GSE74596)[["tx"]]
x <- t(assays(GSE74596_tx)[["TPM"]])
dim(x)
GSE74596$description
GSE74596_tx
table(GSE74596@colData@listData$characteristics_ch1.5)

y_prep<- sub('vα14 inkt thymocyte subset: ', '\\1', GSE74596@colData@listData$characteristics_ch1.5)
y_prep




################
### binomial ###
################
binom <- c(which(y_prep=="NKT0" | y_prep=="NKT17"))
length(binom)
x_mid <- x[binom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #6748
#freq = apply(x_mid, 2, function(x)mean(round(x) > 1))
#mean(freq > 0.25)

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #6,748
y_prep
y <- as.factor(c(y_prep[which(y_prep=="NKT0" | y_prep=="NKT17")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="NKT0" | y_prep=="NKT17"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/GSE74596_NKT0_NKT17.rds")



binom <- c(which(y_prep=="NKT1" | y_prep=="NKT2"))
length(binom)
x_mid <- x[binom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #7879
#freq = apply(x_mid, 2, function(x)mean(round(x) > 1))
#mean(freq > 0.25)

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #6,748
y_prep
y <- as.factor(c(y_prep[which(y_prep=="NKT1" | y_prep=="NKT2")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="NKT1" | y_prep=="NKT2"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/GSE74596_NKT1_NKT2.rds")



binom <- c(which(y_prep=="NKT0" | y_prep=="NKT1"))
length(binom)
x_mid <- x[binom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #5594
#freq = apply(x_mid, 2, function(x)mean(round(x) > 1))
#mean(freq > 0.25)

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #6,748
y_prep
y <- as.factor(c(y_prep[which(y_prep=="NKT0" | y_prep=="NKT1")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="NKT0" | y_prep=="NKT1"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/GSE74596_NKT0_NKT1.rds")



###################
### Multinomial ###
###################


multinom <- c(which(y_prep=="NKT0" | y_prep=="NKT1"| y_prep=="NKT2" | y_prep=="NKT17"))
length(multinom)
x_mid <- x[multinom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #7329 

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) 
y <- as.factor(c(y_prep[which(y_prep=="NKT0" | y_prep=="NKT1"| y_prep=="NKT2" | y_prep=="NKT17")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="NKT0" | y_prep=="NKT1"| y_prep=="NKT2" | y_prep=="NKT17"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/GSE74596_multi.rds")



