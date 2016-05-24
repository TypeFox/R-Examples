### R code from vignette source 'CFR_vignette.Rnw'

###################################################
### code chunk number 1: sourcing
###################################################
library(coarseDataTools)
data(simulated.outbreak.deaths)


###################################################
### code chunk number 2: datapeek
###################################################
simulated.outbreak.deaths[15:20,]


###################################################
### code chunk number 3: preprocessing
###################################################
## set minimum number of observed cases for inclusion
min.cases <- 10

## observed cases
N.1 <- simulated.outbreak.deaths[1:60,"N"]
N.2 <- simulated.outbreak.deaths[61:120,"N"]

## subset to run analyis on times with greater than min.cases
first.t <-  min(which(N.1 > min.cases & N.2 > min.cases))
last.t <-  max(which(N.1 > min.cases & N.2 > min.cases))
idx.for.Estep <- first.t:last.t

## find and label the subset of times to be used for estimation routine
new.times <- 1:length(idx.for.Estep)
simulated.outbreak.deaths <- cbind(simulated.outbreak.deaths, new.times=NA)
simulated.outbreak.deaths[c(idx.for.Estep, idx.for.Estep+60),"new.times"] <- rep(new.times, 2)


###################################################
### code chunk number 4: datapeek2
###################################################
simulated.outbreak.deaths[15:20,]


###################################################
### code chunk number 5: setValues
###################################################
assumed.nu = c(0, .3, .4, .3)
alpha.start <- rep(0, 22)


###################################################
### code chunk number 6: runAnalysis
###################################################
cfr.ests <- EMforCFR(assumed.nu=assumed.nu,
alpha.start.values=alpha.start, full.data=simulated.outbreak.deaths, verb=FALSE,
SEM.var=TRUE, max.iter=500, tol=1e-5)


###################################################
### code chunk number 7: estimationResults
###################################################
cfr.ests$naive.rel.cfr
cfr.ests$glm.rel.cfr
cfr.ests$EM.rel.cfr
cfr.ests$EM.rel.cfr.var.SEM


