### R code from vignette source 'polySegratioMM-overview.Rnw'

###################################################
### code chunk number 1: polySegratioMM-overview.Rnw:63-64
###################################################
library(polySegratioMM)


###################################################
### code chunk number 2: polySegratioMM-overview.Rnw:67-72
###################################################
##library(cacheSweave)  #  comment this out after development phase
## to use:    Sweave("polySegratioMM-overview.Rnw", driver=cacheSweaveDriver)
## or without cache: Sweave("polySegratioMM-overview.Rnw")
op <- options()
options(width=70, digits=4)


###################################################
### code chunk number 3: polySegratioMM-overview.Rnw:90-91
###################################################
data(hexmarkers)


###################################################
### code chunk number 4: polySegratioMM-overview.Rnw:93-100
###################################################
##<<simData, cache=true>>=
## simulate small autohexaploid data set of 500 markers for 200 individuals
## with %70 Single, %20 Double and %10 Triple Dose markers
## created with 
## hexmarkers <- sim.autoMarkers(6,c(0.7,0.2,0.1),n.markers=500,n.individuals=200)
## save(hexmarkers, file="../../data/hexmarkers.RData")
print(hexmarkers)


###################################################
### code chunk number 5: polySegratioMM-overview.Rnw:108-109
###################################################
sr <-  segregationRatios(hexmarkers$markers)


###################################################
### code chunk number 6: polySegratioMM-overview.Rnw:115-118
###################################################
print(plotTheoretical(ploidy.level=6, seg.ratios=sr, main="",
  	    expected.segratio=NULL, proportions=c(0.7,0.2,0.1),
  	    n.individuals=200))


###################################################
### code chunk number 7: polySegratioMM-overview.Rnw:141-148
###################################################
##<<simDataOver, cache=true>>=
## simulate small autohexaploid data set of 500 markers for 200 individuals
## with %70 Single, %20 Double and %10 Triple Dose markers 
## hexmarkers.overdisp <- sim.autoMarkers(6,c(0.7,0.2,0.1),overdispersion=TRUE, shape1=30,n.markers=500,n.individuals=200)
## save(hexmarkers.overdisp, file="../../data/hexmarkers.overdisp.RData")
data(hexmarkers.overdisp)
##print(hexmarkers.overdisp)


###################################################
### code chunk number 8: polySegratioMM-overview.Rnw:154-155
###################################################
sr.overdisp <-  segregationRatios(hexmarkers.overdisp$markers)


###################################################
### code chunk number 9: polySegratioMM-overview.Rnw:160-163
###################################################
print(plotTheoretical(ploidy.level=6, seg.ratios=sr.overdisp, main="",
  expected.segratio=NULL, proportions=c(0.7,0.2,0.1),
  n.individuals=200))


###################################################
### code chunk number 10: polySegratioMM-overview.Rnw:267-268
###################################################
x.mod1 <- setModel(3,6)  # autohexaploid model with 3 components


###################################################
### code chunk number 11: polySegratioMM-overview.Rnw:292-297
###################################################
## produced using the following but loaded as data to avoid the run time on slow machines
##mcmcHexRun <- runSegratioMM(sr.overdisp, x.mod1, burn.in=200, sample=500, plots=FALSE)
##mcmcHexRun <- runSegratioMM(sr.overdisp, x.mod1, plots=FALSE)
## save(mcmcHexRun, file="../../data/mcmcHexRun.RData")
data(mcmcHexRun)


###################################################
### code chunk number 12: polySegratioMM-overview.Rnw:308-309
###################################################
print(mcmcHexRun$run.jags)


###################################################
### code chunk number 13: polySegratioMM-overview.Rnw:313-314
###################################################
print(mcmcHexRun$summary)


###################################################
### code chunk number 14: polySegratioMM-overview.Rnw:321-322
###################################################
print(mcmcHexRun$diagnostics)


###################################################
### code chunk number 15: polySegratioMM-overview.Rnw:326-327
###################################################
print(mcmcHexRun$doses)


###################################################
### code chunk number 16: polySegratioMM-overview.Rnw:338-339
###################################################
print(plot(mcmcHexRun$mcmc.mixture$mcmc.list[[1]][,c("P[1]","mu[1]","sigma","T[140]")]))


###################################################
### code chunk number 17: polySegratioMM-overview.Rnw:364-365
###################################################
print(plot(mcmcHexRun, theoretical=TRUE, main=""))


###################################################
### code chunk number 18: polySegratioMM-overview.Rnw:396-400
###################################################
cat("Employing maximum posterior probability\n")
table(Dose=mcmcHexRun$doses$max.post.dosage, exclude=NULL)
cat("Employing posterior probability > 0.8\n")
table(Dose=mcmcHexRun$doses$dosage[,"0.8"], exclude=NULL)


###################################################
### code chunk number 19: polySegratioMM-overview.Rnw:409-418
###################################################
cat("Employing theChi squared test\n")
dose.chi <- test.segRatio(sr.overdisp, ploidy.level = 6)
table(Chi2Dose=dose.chi$dosage, True=hexmarkers.overdisp$true.doses$dosage, exclude=NULL)
cat("Employing maximum posterior probability\n")
table(MixtureDose=mcmcHexRun$doses$max.post.dosage, True=hexmarkers.overdisp$true.doses$dosage,
exclude=NULL)
cat("Employing posterior probability > 0.8\n")
table(MixtureDose=mcmcHexRun$doses$dosage[,"0.8"], True=hexmarkers.overdisp$true.doses$dosage,
exclude=NULL)


