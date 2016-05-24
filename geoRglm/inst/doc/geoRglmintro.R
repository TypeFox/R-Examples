### R code from vignette source 'geoRglmintro.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: "R settings"
###################################################
options("width"=70)
options(SweaveHooks=list(fig=function() par(mar=c(3,3,1,0.5), mgp=c(2,1,0))))
if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) set.seed(1234)


###################################################
### code chunk number 2: "Loading package" (eval = FALSE)
###################################################
## library(geoR)
## library(geoRglm)


###################################################
### code chunk number 3: geoRglmintro.Rnw:90-94
###################################################
options(geoR.messages = FALSE)
library(geoR)
options(geoR.messages = TRUE)
library(geoRglm)


###################################################
### code chunk number 4: "data-sets"
###################################################
data(b50)
data(p50)


###################################################
### code chunk number 5: geoRglmintro.Rnw:124-125
###################################################
options(geoR.messages = FALSE)


###################################################
### code chunk number 6: geoRglmintro.Rnw:127-129
###################################################
sim.g <- grf(grid = expand.grid(x = seq(1, 10, l = 10), y = seq(1,
10, l = 10)), cov.pars = c(0.1, 0.2))


###################################################
### code chunk number 7: geoRglmintro.Rnw:131-132
###################################################
options(geoR.messages = TRUE)


###################################################
### code chunk number 8: "simulation"
###################################################
sim <- list(coords=sim.g$coords, units.m = c(rep(1, 50), rep(5, 50)))
attr(sim,"class") <- "geodata" 
sim$data <- rpois(100, lambda = sim$units.m*exp(sim.g$data)) 


###################################################
### code chunk number 9: geoRglmintro.Rnw:148-150
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(sim$coords[,1], sim$coords[,2], type = "n")  
text(sim$coords[,1], sim$coords[,2], format(sim$data))                   


###################################################
### code chunk number 10: "Poisson MCMC 1"
###################################################
model2 <- list(cov.pars=c(1,1), beta=1, family="poisson") 
mcmc2.test <- mcmc.control(S.scale=0.2, thin=1)
test2.tune <- glsm.mcmc(p50, model = model2, mcmc.input=mcmc2.test)


###################################################
### code chunk number 11: geoRglmintro.Rnw:192-194
###################################################
mcmc2.tune <- mcmc.control(S.scale=0.5, thin=1)
test2.tune <- glsm.mcmc(p50, model=model2, mcmc.input=mcmc2.tune)


###################################################
### code chunk number 12: geoRglmintro.Rnw:202-204
###################################################
library(coda)
test2.tune.c <- create.mcmc.coda(test2.tune, mcmc.input=mcmc2.tune)


###################################################
### code chunk number 13: geoRglmintro.Rnw:217-221
###################################################
getOption("SweaveHooks")[["fig"]]()
test2.tune.c <- create.mcmc.coda(test2.tune$simulations[45,], mcmc.input=list(S.scale=0.5, thin=1))
par(mfrow=c(1,2))
plot(test2.tune.c, density=FALSE, ask=FALSE, auto.layout = FALSE)
autocorr.plot(test2.tune.c, ask=FALSE, auto.layout = FALSE)


###################################################
### code chunk number 14: "Poisson MCMC 2"
###################################################
mcmc2 <- mcmc.control(S.scale=0.5)
test2 <- glsm.mcmc(p50, model=model2, mcmc.input=mcmc2)


###################################################
### code chunk number 15: "Poisson Kriging 2"
###################################################
out2 <- output.glm.control(sim.predict = TRUE)
pred.test2 <- glsm.krige(test2, locations = cbind(c(0.5,0.5),c(1,0.4)), output = out2)


###################################################
### code chunk number 16: geoRglmintro.Rnw:262-263
###################################################
cbind(pred.test2$predict,pred.test2$mcmc.error)


###################################################
### code chunk number 17: geoRglmintro.Rnw:273-276
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow = c(1,2))
hist(pred.test2$simulations[1,], main="(0.5, 0.5) ")
hist(pred.test2$simulations[2,], main="(1, 0.4)")


###################################################
### code chunk number 18: "Poisson Bayesian Kriging 1"
###################################################
prior5 <- prior.glm.control(phi.prior="fixed", phi=0.1)
mcmc5.tune <- mcmc.control(S.scale=0.01, thin=1)
test5.tune <- pois.krige.bayes(p50, prior=prior5, mcmc.input=mcmc5.tune)


###################################################
### code chunk number 19: "R settings"
###################################################
options(geoR.messages = FALSE)


###################################################
### code chunk number 20: "Poisson Bayesian Kriging 2"
###################################################
mcmc5 <- mcmc.control(S.scale=0.075, thin=100)
out5 <- output.glm.control(threshold=10, quantile=c(0.05,0.99))
test5 <- pois.krige.bayes(p50, locations=t(cbind(c(2.5,3),c(-6050,-3270))), prior=prior5, mcmc.input=mcmc5, output=out5)


###################################################
### code chunk number 21: geoRglmintro.Rnw:349-353
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow = c(1,3))
hist(test5$posterior$simulations[10,], main="(9, 0)")
hist(test5$posterior$simulations[23,], main="(2,2)")
hist(test5$posterior$simulations[36,], main="(5,3)")


###################################################
### code chunk number 22: "Poisson Bayesian Kriging 3"
###################################################
mcmc6.tune <- mcmc.control(S.scale=0.075, n.iter=2000, thin=100, phi.scale=0.01)
prior6 <- prior.glm.control(phi.prior="uniform", phi.discrete=seq(0.02, 1, 0.02), tausq.rel=0.05)
test6.tune <- pois.krige.bayes(p50, prior=prior6, mcmc.input=mcmc6.tune)


###################################################
### code chunk number 23: "Poisson Bayesian Kriging 4"
###################################################
mcmc6 <- mcmc.control(S.scale=0.075, n.iter=400000, thin=200, burn.in=5000, phi.scale=0.12, phi.start=0.5)
test6 <- pois.krige.bayes(p50, locations=t(cbind(c(2.5,3.5),c(-60,-37))), prior=prior6, mcmc.input=mcmc6)


###################################################
### code chunk number 24: geoRglmintro.Rnw:380-384
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow = c(1,3))
hist(test6$posterior$beta$sample, main="beta")
hist(test6$posterior$sigmasq$sample, main="sigmasq")
hist(test6$posterior$phi$sample, main="phi")  


