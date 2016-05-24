### R code from vignette source 'survJamda.Rnw'

###################################################
### code chunk number 1: survJamda.Rnw:28-33 (eval = FALSE)
###################################################
## data(gse4335)
## 
## data(gse4335pheno)
## 
## iter.crossval(gse4335, gse4335pheno[,6], gse4335pheno[,5])


###################################################
### code chunk number 2: survJamda.Rnw:39-40 (eval = FALSE)
###################################################
## iter.crossval(gse4335, gse4335pheno[,6], gse4335pheno[,5], gn.nb =1, gn.nb.display = 1)


###################################################
### code chunk number 3: survJamda.Rnw:44-45 (eval = FALSE)
###################################################
## iter.crossval(gse4335, gse4335pheno[,6], gse4335pheno[,5], plot.roc = 1, gn.nb =1, gn.nb.display = 1)


###################################################
### code chunk number 4: survJamda.Rnw:56-73 (eval = FALSE)
###################################################
## data(gse4335)
## 
## data(gse4335pheno)
## 
## data(gse1992)
## 
## data(gse1992pheno)
## 
## common.gene = intersect(colnames(gse4335),colnames(gse1992))
## 
## data = rbind(gse4335[,common.gene],gse1992[,common.gene])
## 
## surv = c(gse4335pheno[,6],gse1992pheno[,19])
## 
## censor = c(gse4335pheno[,5],gse1992pheno[,18])
## 
## iter.crossval(data,surv,censor,zscore=1)


###################################################
### code chunk number 5: survJamda.Rnw:78-99 (eval = FALSE)
###################################################
## data(gse4335)
## 
## data(gse4335pheno)
## 
## data(gse1992)
## 
## data(gse1992pheno)
## 
## common.gene = intersect(colnames(gse4335),colnames(gse1992))
## 
## data = rbind(gse4335[,common.gene],gse1992[,common.gene])
## 
## surv = c(gse4335pheno[,6],gse1992pheno[,19])
## 
## censor = c(gse4335pheno[,5],gse1992pheno[,18])
## 
## batchID = rep(1,nrow(gse4335))
## 
## batchID = c(batchID,rep(2,nrow(gse1992)))
## 
## iter.crossval.combat(data,surv,censor, batchID)


###################################################
### code chunk number 6: survJamda.Rnw:105-115 (eval = FALSE)
###################################################
## data(gse4335)
## data(gse3143)
## data(gse1992)
## data(gse4335pheno)
## data(gse3143pheno)
## data(gse1992pheno)
## geno.files = c("gse4335","gse3143","gse1992")
## surv.data = list(c(gse4335pheno[,6],gse3143pheno[,4],gse1992pheno[,19]),
## c(gse4335pheno[,5],gse3143pheno[,3],gse1992pheno[,18]))
## main.single.indep.valid(geno.files, surv.data)


###################################################
### code chunk number 7: survJamda.Rnw:119-129 (eval = FALSE)
###################################################
## data(gse4335)
## data(gse3143)
## data(gse1992)
## data(gse4335pheno)
## data(gse3143pheno)
## data(gse1992pheno)
## geno.files = c("gse4335","gse3143","gse1992")
## surv.data = list(c(gse4335pheno[,6],gse3143pheno[,4],gse1992pheno[,19]),
## c(gse4335pheno[,5],gse3143pheno[,3],gse1992pheno[,18]))
## main.merge.indep.valid(geno.files,surv.data)


###################################################
### code chunk number 8: survJamda.Rnw:132-133 (eval = FALSE)
###################################################
## pred.time.indep.valid(geno.files, surv.data)


###################################################
### code chunk number 9: survJamda.Rnw:142-152 (eval = FALSE)
###################################################
## data(gse4335)
## data(gse3143)
## data(gse1992)
## data(gse4335pheno)
## data(gse3143pheno)
## data(gse1992pheno)
## geno.files = c("gse4335","gse3143","gse1992")
## surv.data = list(c(gse4335pheno[,6],gse3143pheno[,4],gse1992pheno[,19]),
## c(gse4335pheno[,5],gse3143pheno[,3],gse1992pheno[,18]))
## meta.main(geno.files, surv.data)


###################################################
### code chunk number 10: survJamda.Rnw:156-157 (eval = FALSE)
###################################################
## proc.simulate()


###################################################
### code chunk number 11: survJamda.Rnw:170-173 (eval = FALSE)
###################################################
## data(gse4335)
## data(gse4335pheno)
## iter.subset(gse4335, gse4335pheno[,6],gse4335pheno[,5])


