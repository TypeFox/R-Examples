### R code from vignette source 'phangorn-specials.Rnw'

###################################################
### code chunk number 1: phangorn-specials.Rnw:48-51
###################################################
options(width=70)
options("show.signif.stars" = FALSE)
foo <- packageDescription("phangorn")


###################################################
### code chunk number 2: phangorn-specials.Rnw:73-79
###################################################
library(phangorn)
data = matrix(c("r","a","y","g","g","a","c","-","c","t","c","g", 
    "a","a","t","g","g","a","t","-","c","t","c","a",                                          
    "a","a","t","-","g","a","c","c","c","t","?","g"), 
    dimnames = list(c("t1", "t2", "t3"),NULL), nrow=3, byrow=TRUE)
data


###################################################
### code chunk number 3: phangorn-specials.Rnw:82-84
###################################################
gapsdata1 = phyDat(data)
gapsdata1


###################################################
### code chunk number 4: phangorn-specials.Rnw:87-90
###################################################
gapsdata2 = phyDat(data, type="USER", levels=c("a","c","g","t","-"), 
    ambiguity = c("?", "n"))
gapsdata2


###################################################
### code chunk number 5: phangorn-specials.Rnw:94-109
###################################################
contrast = matrix(data = c(1,0,0,0,0,
    0,1,0,0,0,
    0,0,1,0,0,
    0,0,0,1,0,   
    1,0,1,0,0,
    0,1,0,1,0,
    0,0,0,0,1,
    1,1,1,1,0,
    1,1,1,1,1),
    ncol = 5, byrow = TRUE)
dimnames(contrast) = list(c("a","c","g","t","r","y","-","n","?"), 
    c("a", "c", "g", "t", "-"))
contrast
gapsdata3 = phyDat(data, type="USER", contrast=contrast)
gapsdata3 


###################################################
### code chunk number 6: phangorn-specials.Rnw:165-170
###################################################
tree = unroot(rtree(3))
fit = pml(tree, gapsdata3)
fit = optim.pml(fit, optQ=TRUE, subs=c(1,0,1,2,1,0,2,1,2,2), 
    control=pml.control(trace=0))
fit


###################################################
### code chunk number 7: phangorn-specials.Rnw:233-243
###################################################
library(phangorn)
primates = read.phyDat("primates.dna", format="phylip", type="DNA")
tree <- NJ(dist.ml(primates))
dat <- phyDat(as.character(primates), "CODON")
fit <- pml(tree, dat)
fit0 <- optim.pml(fit, control = pml.control(trace = 0))
fit1 <- optim.pml(fit, model="codon1", control=pml.control(trace=0))
fit2 <- optim.pml(fit, model="codon2", control=pml.control(trace=0))
fit3 <- optim.pml(fit, model="codon3", control=pml.control(trace=0))
anova(fit0, fit2, fit3, fit1)


###################################################
### code chunk number 8: plotAll
###################################################
trees = allTrees(5)
par(mfrow=c(3,5), mar=rep(0,4)) 
for(i in 1:15)plot(trees[[i]], cex=1, type="u")


###################################################
### code chunk number 9: figAll
###################################################
getOption("SweaveHooks")[["fig"]]()
trees = allTrees(5)
par(mfrow=c(3,5), mar=rep(0,4)) 
for(i in 1:15)plot(trees[[i]], cex=1, type="u")


###################################################
### code chunk number 10: phangorn-specials.Rnw:270-271
###################################################
trees = nni(trees[[1]])


###################################################
### code chunk number 11: phangorn-specials.Rnw:282-283
###################################################
toLatex(sessionInfo())


