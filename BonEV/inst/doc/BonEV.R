### R code from vignette source 'BonEV.Rnw'

###################################################
### code chunk number 1: BonEV.Rnw:40-41 (eval = FALSE)
###################################################
## install.packages("BonEV")


###################################################
### code chunk number 2: BonEV.Rnw:47-48 (eval = FALSE)
###################################################
## library(BonEV)


###################################################
### code chunk number 3: BonEV.Rnw:59-75
###################################################
library(qvalue)
data(hedenfalk)
library(BonEV)
pvalues <- hedenfalk$p
adjp <- Bon_EV(pvalues, 0.05)
summary(adjp)
results <- cbind(adjp$raw_P_value, adjp$BH_adjp, adjp$Storey_adjp, adjp$Bon_EV_adjp)
colnames(results) <- c("raw_P_value", "BH_adjp", "Storey_adjp", "Bon_EV_adjp")
results[1:20,]
summary(results)

##Compare with Benjami-Hochberg and Storey's q-value procedures
sum(adjp$raw_P_value <= 0.05)
sum(adjp$BH_adjp <= 0.05)
sum(adjp$Storey_adjp <= 0.05)
sum(adjp$Bon_EV_adjp <= 0.05)


