
##########################
### Data Set EMTAB2805 ###
##########################

EMTAB2805 <- readRDS("dataset/EMTAB2805.rds")
(EMTAB2805_gene <- experiments(EMTAB2805)[["tx"]]) # tx is possible instead of gene
x <- t(assays(EMTAB2805_gene)[["TPM"]]) # TPM instead of count is possible

y_prep <- as.factor(EMTAB2805@colData@listData$cell_cycle_stage)
y_prep



################
### binomial ###
################
### G1 and G2M ###
binom <- c(which(y_prep=="G1" | y_prep=="G2M"))
length(binom)
x_mid <- x[binom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #13109
#freq = apply(x_mid, 2, function(x)mean(round(x) > 1))
#mean(freq > 0.25)

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #13109
y_prep
y <- as.factor(c(y_prep[which(y_prep=="G1" | y_prep=="G2M")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="G1" | y_prep=="G2M"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/EMTAB2805_G1_G2M.rds")



### S and G1 ###
binom <- c(which(y_prep=="S" | y_prep=="G1"))
length(binom)
x_mid <- x[binom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #12,146
#freq = apply(x_mid, 2, function(x)mean(round(x) > 1))
#mean(freq > 0.25)

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #13109
y_prep
y <- as.factor(c(y_prep[which(y_prep=="S" | y_prep=="G1")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="S" | y_prep=="G1"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/EMTAB2805_S_G1.rds")

### S and G2M ###
binom <- c(which(y_prep=="S" | y_prep=="G2M"))
length(binom)
x_mid <- x[binom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #13,234
#freq = apply(x_mid, 2, function(x)mean(round(x) > 1))
#mean(freq > 0.25)

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) #13109
y_prep
y <- as.factor(c(y_prep[which(y_prep=="S" | y_prep=="G2M")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="S" | y_prep=="G2M"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/EMTAB2805_S_G2M.rds")



###################
### Multinomial ###
###################

multinom <- c(which(y_prep=="S" | y_prep=="G1" | y_prep=="G2M"))
length(multinom)
x_mid <- x[multinom,]
frac_TPM_mid <- rep(NA, dim(x_mid)[2])
for(g in 1:dim(x_mid)[2]){
  frac_TPM_mid[g] <-(sum(round(x_mid[,g]) >1))/dim(x_mid)[1]
}
mean((frac_TPM_mid)>0.25)
table((frac_TPM_mid)>0.25) #12849 

sel_genes <- c(which(frac_TPM_mid >0.25))
length(unique(sel_genes)) 
y <- as.factor(c(y_prep[which(y_prep=="S" | y_prep=="G1" | y_prep=="G2M")]))
y
levels(y) <- 1:(length(levels(y)))
x_sel <- x[which(y_prep=="S" | y_prep=="G1" | y_prep=="G2M"), unique(sel_genes)]
dim(x_sel)
x_sel <- apply(x_sel, 2, as.numeric)
data_aim <- data.frame(x_sel,y)
dim(data_aim)
saveRDS(data_aim, "inst_scRNA/EMTAB2805_multi.rds")






