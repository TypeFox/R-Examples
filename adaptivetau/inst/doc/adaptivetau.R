### R code from vignette source 'adaptivetau.Rnw'

###################################################
### code chunk number 1: adaptivetau.Rnw:8-9
###################################################
pdf.options(pointsize=11) #match latex pointsize


###################################################
### code chunk number 2: adaptivetau.Rnw:68-69
###################################################
library(adaptivetau)


###################################################
### code chunk number 3: adaptivetau.Rnw:72-75
###################################################
transitions = list(c(prey = +1),            # trans 1: prey grows
                   c(prey = -2, pred = +1), # trans 2: predation
                   c(pred = -1))            # trans 3: predator dies


###################################################
### code chunk number 4: adaptivetau.Rnw:78-84
###################################################
lvRateF <- function(x, params, t) {
  return(c(params$r * x["prey"],               # rate of prey growing
           params$beta * x["prey"]*x["pred"] * # rate of predation
             (x["prey"] >= 2),
           params$delta * x["pred"]))          # rate of predators dying
}


###################################################
### code chunk number 5: adaptivetau.Rnw:87-92
###################################################
set.seed(4) # set random number generator seed to be reproducible
simResults = ssa.adaptivetau(init.values = c(prey = 1000, pred = 500),
                             transitions, lvRateF,
                             params = list(r=10, beta=0.01, delta=10),
                             tf=12)


###################################################
### code chunk number 6: adaptivetau.Rnw:95-98 (eval = FALSE)
###################################################
## matplot(simResults[,"time"], simResults[,c("prey","pred")], type='l',
##         xlab='Time', ylab='Counts (log scale)', log='y')
## legend("bottomleft", legend=c("prey", "predator"), lty=1:2, col=1:2)


###################################################
### code chunk number 7: adaptivetau.Rnw:103-107
###################################################
par(mex=1, mar=c(3,3,1,1), mgp=c(2,.75,0))
matplot(simResults[,"time"], simResults[,c("prey","pred")], type='l',
        xlab='Time', ylab='Counts (log scale)', log='y')
legend("bottomleft", legend=c("prey", "predator"), lty=1:2, col=1:2)


###################################################
### code chunk number 8: adaptivetau.Rnw:114-118 (eval = FALSE)
###################################################
## simResults = ssa.exact(init.values = c(prey = 1000, pred = 500),
##                        transitions, lvRateF,
##                        params = list(r=10, beta=0.01, delta=10),
##                        tf=12)


###################################################
### code chunk number 9: adaptivetau.Rnw:121-126 (eval = FALSE)
###################################################
## simResults = ssa.adaptivetau(init.values = c(prey = 1000, pred = 500),
##                              transitions, lvRateF,
##                              params = list(r=10, beta=0.01, delta=10),
##                              tf=12,
##                              tl.params = list(epsilon = 0.005))


###################################################
### code chunk number 10: adaptivetau.Rnw:139-147
###################################################
library(adaptivetau)
transitions = cbind(c(-1,+1), # num individuals w/ derived allele decreases
                    c(+1,-1)) # num individuals w/ derived allele increases
driftRateF <- function(x, params, t) {
  rep(x[1] * x[2]/sum(x)^2, 2)
}
set.seed(1) # set random number generator seed to be reproducible
r=ssa.adaptivetau(c(500,500), transitions, driftRateF, params=NULL, tf=Inf)


###################################################
### code chunk number 11: adaptivetau.Rnw:152-154
###################################################
par(mex=1, mar=c(3,3,1,1), mgp=c(2,.75,0))
plot(r[,"time"], r[,"x1"], type='l', xlab='Time', ylab='Number of derived alleles', ylim=c(0,1000))


###################################################
### code chunk number 12: adaptivetau.Rnw:175-205
###################################################
library(adaptivetau)
init.values = c(
  S = 10^5, # susceptible humans
  I1 = 0,   # infected humans
  I2 = 0,   # infected humans
  R = 0)    # recovered (and immune) humans

transitions = list(
    c(S = -1, I1 = +1), # infection (animal-adapted strain)
    c(S = -1, I2 = +1), # infection (human-adapted strain)
    c(I1 = -1),         # death due to infection
    c(I2 = -1),
    c(I1 = -1, R = +1), # recovery
    c(I2 = -1, R = +1)
    )

SIRrateF <- function(x, p, t) {
  return(c(x["S"] * (p$zoonotic + p$beta[1]*x["I1"]), # infection rate
           x["S"] * (p$beta[2]*x["I2"] + p$mu*p$beta[1]*x["I1"]),
           params$delta[1]*x["I1"], # infected death rate
           params$delta[2]*x["I2"],
           params$gamma[1]*x["I1"], # recovery rate
           params$gamma[2]*x["I2"]))
}

params = list(zoonotic=1e-6, beta=c(1e-7, 1e-5), mu=0.1,
              delta=c(1e-2,1e-5), gamma=c(0.1, 0.1))

set.seed(3) # set random number generator seed to be reproducible
r=ssa.adaptivetau(init.values, transitions, SIRrateF, params, tf=1000)


###################################################
### code chunk number 13: adaptivetau.Rnw:210-213
###################################################
par(mex=1, mar=c(3,3,1,1), mgp=c(2,.75,0))
matplot(r[,"time"], r[,c("S","I1","I2")], type='l', log='y', xlab="Time", ylab="Individuals (log scale)", col=c('black', 'red', 'blue'))
legend("topright", legend=c("S", expression(I[1]), expression(I[2])), lty=1:3, col=c('black', 'red', 'blue'))


