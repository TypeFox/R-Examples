### R code from vignette source 'hwde.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(show.signif.stars=FALSE)


###################################################
### code chunk number 2: chunk2
###################################################
obs <- c(AA=147, Aa=78, aa=17)
oset <- c(0, log(2), 0)
ma <- c(0,1,2)
hw.glm <- glm(obs ~ ma, family=poisson, offset=oset)
summary(hw.glm)


###################################################
### code chunk number 3: chunk3
###################################################
hwdat <- data.frame(Observed=c(147,78,17), locus1=c("AA","Aa","aa"))


###################################################
### code chunk number 4: chunk4
###################################################
library(hwde)
hwde(data=hwdat)


###################################################
### code chunk number 5: chunk5-maa
###################################################
data.df <- hwde(data = hwdat)$data.df
names(data.df)
summary(glm(obs ~ a + aa, offset=oset, family=poisson, data=data.df))$coef


###################################################
### code chunk number 6: chunk6 (eval = FALSE)
###################################################
## hwde(data=mendelABC, loci=c("seedshape","cotylcolor","coatcolor"))


###################################################
### code chunk number 7: chunk7 (eval = FALSE)
###################################################
## hwdat <- read.table("hw.txt", header=TRUE)


###################################################
### code chunk number 8: chunk8 (eval = FALSE)
###################################################
## IndianIrish <- read.table("IndianIrish.txt", header=TRUE)


###################################################
### code chunk number 9: chunk9
###################################################
hwde(data=IndianIrish)


###################################################
### code chunk number 10: chunk10 (eval = FALSE)
###################################################
## hwde(data=IndianIrish, group.terms=FALSE)


###################################################
### code chunk number 11: chunk11-all9
###################################################
II.hwde <- hwde(data = mendelABC, loci = c("seedshape", "cotylcolor",
              "coatcolor"), keep.models=T)
models <- II.hwde$models
maxmodel <- models[[length(models)]]
summary(maxmodel)$coef


###################################################
### code chunk number 12: chunk12
###################################################
hwdat.hw <- hwde(data=hwdat)
names(hwdat)
hwdat.hw$data.df


###################################################
### code chunk number 13: chunk13
###################################################
data.df <- hwdat.hw$data.df
model1 <- glm(obs ~ a, family=poisson, data=data.df, offset=log(oset))
model2 <- glm(obs ~ a+aa, family=poisson, data=data.df, offset=log(oset))
model1


###################################################
### code chunk number 14: chunk14
###################################################
II.hw <- hwde(data=IndianIrish, aovtable.print=FALSE)
dataII.df <- II.hw$data.df
dataII.df


###################################################
### code chunk number 15: chunk15
###################################################
hwde(termlist=c("+a","+aa"), refmodel=c(1,2), data=hwdat)


###################################################
### code chunk number 16: chunk16
###################################################
hwde(termlist=c("+gp","+a","+b","+a+b","+aa"), refmodel=c(1,2,2,2,5),
     data=IndianIrish)


###################################################
### code chunk number 17: chunk17
###################################################
hwdat.hw <- hwde(data=hwdat, keep.models=TRUE)
hwdat.hw$models[[2]]          # The Hardy-Weinberg model
fitted(hwdat.hw$models[[2]])


