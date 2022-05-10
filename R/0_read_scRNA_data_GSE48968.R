GSE48968 <- readRDS("dataset/GSE48968-GPL13112.rds")
(GSE48968_gene <- experiments(GSE48968)[["gene"]])
data <- t(assays(GSE48968_gene)[["TPM"]])
dim(data)
y_org <- as.factor(GSE48968@colData@listData$characteristics_ch1.2)
y_org<- sub('stimulation: ', '\\1', y_org)
y_org
table(y_org)

y_prep <- substr(y_org, 1, 6)
y <- as.factor(y_prep)
levels(y) <- 1:(length(levels(y)))
x <- apply(data[, 1:(dim(data)[2])], 2, as.numeric)
dim(x)
data_aim <- data.frame(x,y)
saveRDS(data_aim, "inst_scRNA/GSE48968.rds")



################
### binomial ###
################
### 1h LPS and 4h LPS ###
binom <- c(which(y_prep=="1h LPS" | y_prep=="4h LPS"))
length(binom)
x_mid <- x[binom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #7987 for 1h LPS and 4h LPS
#freq = apply(x_mid, 2, function(x)mean(round(x) > 1))
#mean(freq > 0.25)

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #7987
y <- as.factor(c(y_prep[which(y_prep=="1h LPS" | y_prep=="4h LPS")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="1h LPS" | y_prep=="4h LPS"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/GSE48968_1h_4h.rds")


### 4h LPS and 6h LPS ###
binom <- c(which(y_prep=="4h LPS" | y_prep=="6h LPS"))
length(binom)
x_mid <- x[binom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #7370 for 4h LPS and 6h LPS
#freq = apply(x_mid, 2, function(x)mean(round(x) > 1))
#mean(freq > 0.25)

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #7370
y <- as.factor(c(y_prep[which(y_prep=="4h LPS" | y_prep=="6h LPS")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="4h LPS" | y_prep=="6h LPS"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)

saveRDS(data_aim, "inst_scRNA/GSE48968_4h_6h.rds")



###################
### multinomial ###
###################

multinom <- c(which(y_prep=="1h LPS" | y_prep=="4h LPS" | y_prep=="6h LPS"))
length(multinom)
x_mid <- x[multinom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #7831 

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #11065
y <- as.factor(c(y_prep[which(y_prep=="1h LPS" | y_prep=="4h LPS" | y_prep=="6h LPS")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="1h LPS" | y_prep=="4h LPS" | y_prep=="6h LPS"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)

saveRDS(data_aim, "inst_scRNA/GSE48968_multi.rds")






