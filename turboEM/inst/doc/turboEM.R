### R code from vignette source 'turboEM.Rnw'

###################################################
### code chunk number 1: load
###################################################
library(turboEM)


###################################################
### code chunk number 2: help (eval = FALSE)
###################################################
## help(package="turboEM")


###################################################
### code chunk number 3: data
###################################################
poissmix.dat <- data.frame(death=0:9,
			   freq=c(162,267,271,185,111,61,27,8,3,1))
y <- poissmix.dat$freq


###################################################
### code chunk number 4: fixptfn
###################################################
fixptfn <- function(p, y) {
	pnew <- rep(NA,3)
	i <- 0:(length(y)-1)
	denom <- p[1]*exp(-p[2])*p[2]^i + (1 - p[1])*exp(-p[3])*p[3]^i
	zi <- p[1]*exp(-p[2])*p[2]^i / denom
	pnew[1] <- sum(y*zi)/sum(y)
	pnew[2] <- sum(y*i*zi)/sum(y*zi)
	pnew[3] <- sum(y*i*(1-zi))/sum(y*(1-zi))
	p <- pnew
	return(pnew)
}


###################################################
### code chunk number 5: objfn
###################################################
objfn <- function(p, y) {
	i <- 0:(length(y)-1)
	loglik <- y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) +
			(1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
	return ( -sum(loglik) )
}


###################################################
### code chunk number 6: fit1
###################################################
res <- turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
	       method=c("em", "squarem", "pem"), y=y)
options(digits=13)
res


###################################################
### code chunk number 7: pars
###################################################
pars(res)


###################################################
### code chunk number 8: showmethods
###################################################
options(digits=7)
grad(res)
hessian(res)
stderror(res)


###################################################
### code chunk number 9: plot
###################################################
res1 <- turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
		method=c("em", "squarem", "pem"), y=y,
		control.run=list(keep.objfval=TRUE))
res1
plot(res1, xlim=c(0.001, 0.02))


###################################################
### code chunk number 10: boundary
###################################################
boundary <- function(par, dr) {
	lower <- c(0, 0, 0)
	upper <- c(1, 10000, 10000)
	low1 <- max(pmin((lower-par)/dr, (upper-par)/dr))
	upp1 <- min(pmax((lower-par)/dr, (upper-par)/dr))
	return(c(low1, upp1))
}


###################################################
### code chunk number 11: fit1
###################################################
res2 <- turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
		boundary=boundary, method="decme", y=y)
options(digits=13)
res2


###################################################
### code chunk number 12: noobjfn
###################################################
res3 <- turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, boundary=boundary, y=y)
res3


###################################################
### code chunk number 13: errorcall
###################################################
error(res3)


###################################################
### code chunk number 14: noobjfn
###################################################
res4 <- turboem(par=c(0.9, 1, 3), fixptfn=fixptfn, objfn=objfn,
		boundary=boundary, y=y)
res4


###################################################
### code chunk number 15: err
###################################################
error(res4)


###################################################
### code chunk number 16: fit3
###################################################
pconstr <- function(par) {
	lower <- c(0, 0, 0)
	upper <- c(1, Inf, Inf)
	return(all(lower < par & par < upper))
}
res5 <- turboem(par=c(0.9, 1, 3), fixptfn=fixptfn, objfn=objfn,
		boundary=boundary, y=y, pconstr=pconstr)
res5


###################################################
### code chunk number 17: changetol
###################################################
res6 <- turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
		method=c("em", "pem", "squarem"), y=y,
		control.run=list(tol=1.0e-10))
res6


###################################################
### code chunk number 18: objfnconv
###################################################
res7 <- turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
		method=c("em", "pem", "squarem"), y=y,
		control.run=list(tol=1.0e-10, convtype="objfn"))
res7


###################################################
### code chunk number 19: userdefconv
###################################################
convfn.user <- function(old, new) {
	max(abs(new-old)) < tol
}
res8 <- turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
		method=c("em", "pem", "squarem"), y=y,
		control.run=list(tol=1.0e-10, convfn.user = convfn.user))
res8


###################################################
### code chunk number 20: userdefconv2
###################################################
convfn.user.objfn <- function(old, new) {
	abs(new - old)/(abs(old) + 1) < tol
}
res9 <- turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
		method=c("em", "pem", "squarem"), y=y,
		control.run=list(tol=1.0e-8, convtype="objfn",
		convfn.user = convfn.user.objfn))
res9


###################################################
### code chunk number 21: newstop
###################################################
res10 <- turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
		method=c("em", "pem", "squarem"), y=y,
		control.run=list(tol=1.0e-15, stoptype="maxtime",
		maxtime=10))
res10


###################################################
### code chunk number 22: newstop2
###################################################
stopfn.user <- function() {
	iter >= maxiter | elapsed.time >= maxtime
}
res11 <- turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
		method=c("em", "pem", "squarem"), y=y,
		control.run=list(tol=1.0e-15, stopfn.user=stopfn.user,
		maxtime=0.2, maxiter=2000))
res11


###################################################
### code chunk number 23: fit2
###################################################
res12 <- turboem(par = c(0.9, 1, 3), fixptfn=fixptfn, objfn=objfn,
		 boundary=boundary, pconstr=pconstr,
		 method=c("em", "squarem", "squarem", "decme", "decme",
		          "qn", "qn", "pem", "pem"),
		 control.method=list(list(), list(K=2), list(K=3),
		     list(version=2), list(version="2s"),
		     list(qn=1), list(qn=2),
		     list(version="arithmetic"), list(version="geometric")),
		 y=y)
res12


###################################################
### code chunk number 24: parallelMCregister
###################################################
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)


###################################################
### code chunk number 25: runMC
###################################################
time.parallel <- system.time(res.parallel <-
	turboem(par=c(0.9, 1, 6), fixptfn=fixptfn, objfn=objfn,
		method=c("em", "pem", "squarem"), y=y, parallel=TRUE,
		control.run=list(tol=1.0e-14, stoptype="maxtime",
		maxtime=10)))
res.parallel
time.parallel


###################################################
### code chunk number 26: stopParallel
###################################################
stopCluster(cl)


###################################################
### code chunk number 27: seqCompare
###################################################
time.sequential <- system.time(res.sequential <-
	turboem(par=c(0.9, 1, 6), fixptfn=fixptfn, objfn=objfn,
		method=c("em", "pem", "squarem"), y=y, parallel=FALSE,
		control.run=list(tol=1.0e-14, stoptype="maxtime",
		maxtime=10)))
res.sequential
time.sequential


###################################################
### code chunk number 28: control
###################################################
method.names <- c("EM", "squaremK1", "squaremK2", "parabolicEM",
		  "dynamicECME", "quasiNewton")
nmethods <- length(method.names)
method <- c("em", "squarem", "squarem", "pem", "decme", "qn")
control.method <- vector("list", nmethods)
names(control.method) <- method.names
control.method[["EM"]] <- list()
control.method[["squaremK1"]] <- list(K=1)
control.method[["squaremK2"]] <- list(K=2)
control.method[["parabolicEM"]] <- list(version="geometric")
control.method[["dynamicECME"]] <- list(version="2s")
control.method[["quasiNewton"]] <- list(qn=2)


###################################################
### code chunk number 29: runparams
###################################################
control.run <- list(tol=1e-7, stoptype="maxtime", maxtime=2,
		    convtype="parameter")


###################################################
### code chunk number 30: seed
###################################################
NREP <- 100
library(setRNG)
test.rng <- list(kind = "Mersenne-Twister",
		 normal.kind = "Inversion", seed = 1)
setRNG(test.rng)
starting.values <- cbind(runif(NREP),runif(NREP,0,4),runif(NREP,0,4))
head(starting.values, 3)


###################################################
### code chunk number 31: run
###################################################
simtime <- system.time(
     results <- turboSim(parmat=starting.values, fixptfn=fixptfn,
		    objfn=objfn, method=method, boundary=boundary,
		    pconstr=pconstr, method.names=method.names,
		    y=y, control.method=control.method,
		    control.run=control.run)
		       )
simtime


###################################################
### code chunk number 32: runParallel
###################################################
cl <- makeCluster(2)
simtime.par <- system.time(
     results.par <- turboSim(parmat=starting.values, fixptfn=fixptfn,
		    objfn=objfn, method=method, boundary=boundary,
		    pconstr=pconstr, method.names=method.names,
		    y=y, control.method=control.method,
		    control.run=control.run, parallel=TRUE)
			   )
simtime.par
stopCluster(cl)


###################################################
### code chunk number 33: turboSimPrint
###################################################
class(results)
results


###################################################
### code chunk number 34: table
###################################################
summary(results, eps=0.01)


###################################################
### code chunk number 35: table2
###################################################
summary(results, eps=0.01, sol=1989.945859883)


###################################################
### code chunk number 36: boxplots
###################################################
fails <- with(results, fail | !convergence |
              value.objfn > apply(value.objfn, 1, min) + 0.01)
boxplot(results, whichfail=fails)


###################################################
### code chunk number 37: dataprofile
###################################################
dataprof(results)


###################################################
### code chunk number 38: scatterplot
###################################################
pairs(results, which.methods=1:4, cex=0.8, whichfail=fails)


