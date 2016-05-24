### R code from vignette source 'STERGM.Snw'

###################################################
### code chunk number 1: STERGM.Snw:25-26
###################################################
options(width=60)


###################################################
### code chunk number 2: STERGM.Snw:29-30
###################################################
options(continue=" ")


###################################################
### code chunk number 3: STERGM.Snw:33-34
###################################################
library('tergm')


###################################################
### code chunk number 4: STERGM.Snw:55-59 (eval = FALSE)
###################################################
## install.packages("ergm")
## install.packages("networkDynamic")
## library(ergm)
## library(networkDynamic)


###################################################
### code chunk number 5: STERGM.Snw:89-92
###################################################
library(ergm)
data("florentine")
ls()


###################################################
### code chunk number 6: flobusplot
###################################################
plot(flobusiness)


###################################################
### code chunk number 7: STERGM.Snw:110-112
###################################################
fit1 <- ergm(flobusiness~edges+gwesp(0,fixed=T))
summary(fit1)


###################################################
### code chunk number 8: sim1plot
###################################################
sim1 <- simulate(fit1,nsim=1,
          control=control.simulate.ergm(MCMC.burnin=1000))
plot(sim1)


###################################################
### code chunk number 9: STERGM.Snw:324-325
###################################################
theta.diss <- log(9)


###################################################
### code chunk number 10: STERGM.Snw:334-342
###################################################

stergm.fit.1 <- stergm(flobusiness,
	formation= ~edges+gwesp(0,fixed=T),
	dissolution = ~offset(edges),
	targets="formation",
	offset.coef.diss = theta.diss,
	estimate = "EGMME"
	)


###################################################
### code chunk number 11: fit1diag
###################################################
mcmc.diagnostics(stergm.fit.1)


###################################################
### code chunk number 12: STERGM.Snw:368-373
###################################################
stergm.fit.1
names(stergm.fit.1)
stergm.fit.1$formation
stergm.fit.1$formation.fit
summary(stergm.fit.1)


###################################################
### code chunk number 13: STERGM.Snw:383-385
###################################################
stergm.sim.1 <- simulate.stergm(stergm.fit.1, nsim=1, 
    time.slices = 1000)


###################################################
### code chunk number 14: STERGM.Snw:401-402
###################################################
stergm.sim.1


###################################################
### code chunk number 15: STERGM.Snw:427-428
###################################################
network.extract(stergm.sim.1,at=429)


###################################################
### code chunk number 16: simex
###################################################
plot(network.extract(stergm.sim.1,at=882))


###################################################
### code chunk number 17: STERGM.Snw:447-449
###################################################
summary(flobusiness~edges+gwesp(0,fixed=T))
colMeans(attributes(stergm.sim.1)$stats)


###################################################
### code chunk number 18: statsform1
###################################################
plot(attributes(stergm.sim.1)$stats)


###################################################
### code chunk number 19: STERGM.Snw:469-472
###################################################
stergm.sim.1.dm <- as.data.frame(stergm.sim.1)
names(stergm.sim.1.dm)
mean(stergm.sim.1.dm$duration)


###################################################
### code chunk number 20: STERGM.Snw:483-484
###################################################
get.edge.value(stergm.sim.1, "active", unlist=FALSE)[[25]]


###################################################
### code chunk number 21: STERGM.Snw:543-544
###################################################
theta.diss.100 <- log(99)


###################################################
### code chunk number 22: STERGM.Snw:550-553
###################################################
summary(fit1)
theta.form <- fit1$coef 
theta.form


###################################################
### code chunk number 23: STERGM.Snw:559-560
###################################################
theta.form[1] <- theta.form[1] - theta.diss.100


###################################################
### code chunk number 24: STERGM.Snw:565-572
###################################################
stergm.sim.2 <- simulate(flobusiness,
	formation=~edges+gwesp(0,fixed=T),
	dissolution=~edges,
	monitor="all",
	coef.form=theta.form,
	coef.diss=theta.diss.100,
	time.slices=10000)


###################################################
### code chunk number 25: simform
###################################################
summary(flobusiness~edges+gwesp(0,fixed=T))
colMeans(attributes(stergm.sim.2)$stats)
stergm.sim.dm.2 <- as.data.frame(stergm.sim.2)
mean(stergm.sim.dm.2$duration)
plot(attributes(stergm.sim.2)$stats)


###################################################
### code chunk number 26: STERGM.Snw:601-603
###################################################
data(samplk)
ls(pattern="samp*")


###################################################
### code chunk number 27: STERGM.Snw:607-610
###################################################
samp <- list()
samp[[1]] <- samplk1
samp[[2]] <- samplk2


###################################################
### code chunk number 28: STERGM.Snw:616-617
###################################################
plot(samplk1)


###################################################
### code chunk number 29: STERGM.Snw:632-637
###################################################
stergm.fit.3 <- stergm(samp,
	formation= ~edges+mutual+ctriad+ttriad,
	dissolution = ~edges+mutual+ctriad+ttriad,
	estimate = "CMLE"
	)


###################################################
### code chunk number 30: STERGM.Snw:647-648
###################################################
summary(stergm.fit.3)


###################################################
### code chunk number 31: STERGM.Snw:685-688
###################################################
msm.net <- network.initialize(500, directed=F)	
msm.net %v% 'race' <- c(rep(0,250),rep(1,250))
msm.net


###################################################
### code chunk number 32: STERGM.Snw:694-697
###################################################
msm.form.formula <- ~edges+nodematch('race')+degree(0)+
    concurrent
msm.target.stats <- c(225,187,180,90)


###################################################
### code chunk number 33: STERGM.Snw:707-708
###################################################
msm.diss.formula <- ~offset(edges)+offset(nodematch("race"))


###################################################
### code chunk number 34: STERGM.Snw:736-737
###################################################
msm.theta.diss <- c(2.944, -0.747) 


###################################################
### code chunk number 35: STERGM.Snw:747-756
###################################################
set.seed(0)
msm.fit <- stergm(msm.net,
	formation= msm.form.formula,
	dissolution= msm.diss.formula,
	targets="formation",
	target.stats= msm.target.stats,
	offset.coef.diss = msm.theta.diss,
	estimate = "EGMME"
)


###################################################
### code chunk number 36: msmdiag
###################################################
mcmc.diagnostics(msm.fit)


###################################################
### code chunk number 37: STERGM.Snw:777-778
###################################################
summary(msm.fit)


###################################################
### code chunk number 38: STERGM.Snw:783-784
###################################################
msm.sim <- simulate(msm.fit,time.slices=1000)


###################################################
### code chunk number 39: STERGM.Snw:789-791
###################################################
colMeans(attributes(msm.sim)$stats)
msm.target.stats


###################################################
### code chunk number 40: msmht
###################################################
msm.sim.dm <- as.data.frame(msm.sim)
plot(msm.sim.dm$head,msm.sim.dm$tail)


###################################################
### code chunk number 41: STERGM.Snw:809-817
###################################################
names(msm.sim.dm)
msm.sim.dm$race1 <- msm.sim.dm$head>250
msm.sim.dm$race2 <- msm.sim.dm$tail>250
msm.sim.dm$homoph <- msm.sim.dm$race1 == msm.sim.dm$race2
mean(msm.sim.dm$duration[msm.sim.dm$homoph==T & 
  msm.sim.dm$onset.censored==F & msm.sim.dm$terminus.censored==F ])
mean(msm.sim.dm$duration[msm.sim.dm$homoph==F & 
  msm.sim.dm$onset.censored==F & msm.sim.dm$terminus.censored==F ])


