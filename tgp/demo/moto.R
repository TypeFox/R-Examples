###################################################
### chunk number 1: 
###################################################
library(tgp)
##options(width=65)
seed <- 0; set.seed(seed)


###################################################
### chunk number 2: 
###################################################
library(MASS)
X <- data.frame(times=mcycle[,1])
Z <- data.frame(accel=mcycle[,2])


###################################################
### chunk number 3: 
###################################################
moto.bgp <- bgp(X=X, Z=Z, verb=0)


###################################################
### chunk number 4: bgp
###################################################
plot(moto.bgp, main='GP,', layout='surf')


###################################################
### chunk number 5: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


###################################################
### chunk number 6: 
###################################################
moto.btlm <- btlm(X=X, Z=Z, verb=0)


###################################################
### chunk number 7: btlm
###################################################
plot(moto.btlm, main='Bayesian CART,', layout='surf')


###################################################
### chunk number 8: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


###################################################
### chunk number 9: 
###################################################
moto.btgpllm <- btgpllm(X=X, Z=Z, bprior="b0", verb=0)
moto.btgpllm.p <- predict(moto.btgpllm) ## using MAP


###################################################
### chunk number 10: btgp
###################################################
par(mfrow=c(1,2))
plot(moto.btgpllm, main='treed GP LLM,', layout='surf')
plot(moto.btgpllm.p, center='km', layout='surf')


###################################################
### chunk number 11: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


###################################################
### chunk number 12: btgpq
###################################################
par(mfrow=c(1,2))
plot(moto.btgpllm, main='treed GP LLM,', layout='as')
plot(moto.btgpllm.p, as='ks2', layout='as')


###################################################
### chunk number 13: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


