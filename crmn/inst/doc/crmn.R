### R code from vignette source 'crmn.Rnw'

###################################################
### code chunk number 1: load
###################################################
library(crmn)
help(package="crmn")


###################################################
### code chunk number 2: mix-example
###################################################
data(mix)
head(fData(mix))[,1:4]
head(pData(mix))


###################################################
### code chunk number 3: div
###################################################
Ys <- standards(mix)
Ya <- analytes(mix)
dim(Ys)
dim(Ya)


###################################################
### code chunk number 4: nfit-eset
###################################################
nfit <- normFit(mix, "crmn", factor="type", ncomp=2)


###################################################
### code chunk number 5: nfit-matrix
###################################################
G <- model.matrix(~-1+mix$type)
nfit <- normFit(mix, "crmn", factor=G, ncomp=2)


###################################################
### code chunk number 6: nfit-eset-w-q2
###################################################
nfit <- normFit(mix, "crmn", factor="type")
#complexty (number of PC's):
sFit(nfit)$ncomp


###################################################
### code chunk number 7: tz
###################################################
slplot(sFit(nfit)$fit$pc, scol=as.integer(mix$runorder))


###################################################
### code chunk number 8: plot
###################################################
nfit
plot(nfit)


###################################################
### code chunk number 9: norm
###################################################
normed.crmn <- normPred(nfit, mix, factor="type")


###################################################
### code chunk number 10: alternative
###################################################
normed.one <- normalize(mix, "one", one="Hexadecanoate_13C4")
normed.nomis <- normalize(mix, "nomis")


###################################################
### code chunk number 11: compare
###################################################

pca.crmn <- pca(scale(log(t(exprs(normed.crmn)))))
pca.one <- pca(scale(log(t(exprs(normed.one)))))
pca.nomis <- pca(scale(log(t(exprs(normed.nomis)))))
par(mfrow=c(1,3))
plot(scores(pca.one), col=as.integer(mix$type),
     pch=as.integer(mix$runorder),
     main="Single IS")
plot(scores(pca.nomis), col=as.integer(mix$type),
     pch=as.integer(mix$runorder),
     main="NOMIS")
plot(scores(pca.crmn), col=as.integer(mix$type),
     pch=as.integer(mix$runorder),
     main="CRMN")



###################################################
### code chunk number 12: mix-example
###################################################
Y <- exprs(mix)
replicates <- factor(mix$type)
G <- model.matrix(~-1+replicates)
isIS <- fData(mix)$tag == 'IS'


###################################################
### code chunk number 13: div
###################################################
standards(Y, isIS)
analytes(Y, isIS)


###################################################
### code chunk number 14: nfit1
###################################################
nfit <- normFit(Y, "crmn", factors=G, ncomp=2, standards=isIS)


###################################################
### code chunk number 15: norm
###################################################
normed.crmn <- normPred(nfit, Y, factors=G, standards=isIS, ncomp=2)


###################################################
### code chunk number 16: alternative
###################################################
normed.crmn <- normalize(Y, "crmn", factors=G, standards=isIS, ncomp=2)


