### R code from vignette source 'sim_condBD.Rnw'

###################################################
### code chunk number 1: preamble
###################################################
options(continue="+");


###################################################
### code chunk number 2: sourceCode
###################################################
library(DOBAD)


###################################################
### code chunk number 3: sim_condBD.Rnw:170-175
###################################################
L <- .3; m <- .5; 
nu <- .4
set.seed(112)
unobservedChain <- birth.death.simulant(t=5, X0=11, lambda=.3, mu=.5, nu=.4);
unobservedChain;


###################################################
### code chunk number 4: sim_condBD.Rnw:178-181
###################################################
times <- c(0, .21,.62,.73, 1.44, 1.95, 3.56, 4.17);
obsData <- getPartialData(times,  unobservedChain);
obsData;


###################################################
### code chunk number 5: condSim
###################################################
nsims <- 2;
condSims <- sim.condBD(N=nsims, bd.PO=obsData, L=L, m=m, nu=nu);
condSims[1]
condSims[2]


###################################################
### code chunk number 6: conclusion
###################################################
options(continue=" "); ##undo what we set at top


