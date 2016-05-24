### R code from vignette source 'Ancestral.Rnw'

###################################################
### code chunk number 1: Ancestral.Rnw:44-47
###################################################
options(width=70)
options("show.signif.stars" = FALSE)
foo <- packageDescription("phangorn")


###################################################
### code chunk number 2: Ancestral.Rnw:62-67
###################################################
library(phangorn)
primates = read.phyDat("primates.dna", format = "phylip", type = "DNA")
tree = pratchet(primates, trace=0)
tree = acctran(tree, primates) 
parsimony(tree, primates)


###################################################
### code chunk number 3: Ancestral.Rnw:73-75
###################################################
anc.acctran = ancestral.pars(tree, primates, "ACCTRAN")
anc.mpr = ancestral.pars(tree, primates, "MPR")


###################################################
### code chunk number 4: plotLOGO
###################################################
tmp <- require(seqLogo)
if(tmp) seqLogo( t(subset(anc.mpr, getRoot(tree), 1:20)[[1]]), ic.scale=FALSE)


###################################################
### code chunk number 5: figLOGO
###################################################
getOption("SweaveHooks")[["fig"]]()
tmp <- require(seqLogo)
if(tmp) seqLogo( t(subset(anc.mpr, getRoot(tree), 1:20)[[1]]), ic.scale=FALSE)


###################################################
### code chunk number 6: Ancestral.Rnw:94-96
###################################################
options(SweaveHooks=list(fig=function()
par(mar=c(2.1, 4.1, 2.1, 2.1))))


###################################################
### code chunk number 7: plotMPR
###################################################
par(mfrow=c(2,1))
plotAnc(tree, anc.mpr, 17)
title("MPR")
plotAnc(tree, anc.acctran, 17)
title("ACCTRAN")


###################################################
### code chunk number 8: figMPR
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,1))
plotAnc(tree, anc.mpr, 17)
title("MPR")
plotAnc(tree, anc.acctran, 17)
title("ACCTRAN")


###################################################
### code chunk number 9: Ancestral.Rnw:123-125
###################################################
fit = pml(tree, primates)
fit = optim.pml(fit, model="F81", control = pml.control(trace=0))


###################################################
### code chunk number 10: Ancestral.Rnw:137-139
###################################################
anc.ml = ancestral.pml(fit, "ml")
anc.bayes = ancestral.pml(fit, "bayes")


###################################################
### code chunk number 11: plotMLB
###################################################
par(mfrow=c(2,1))
plotAnc(tree, anc.ml, 17)
title("ML")
plotAnc(tree, anc.bayes, 17)
title("Bayes")


###################################################
### code chunk number 12: figMLB
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,1))
plotAnc(tree, anc.ml, 17)
title("ML")
plotAnc(tree, anc.bayes, 17)
title("Bayes")


###################################################
### code chunk number 13: Ancestral.Rnw:163-164
###################################################
toLatex(sessionInfo())


