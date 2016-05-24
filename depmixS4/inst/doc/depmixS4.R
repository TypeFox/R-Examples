### R code from vignette source 'depmixS4.Rnw'

###################################################
### code chunk number 1: depmixS4.Rnw:72-74
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE, digits = 4)
library("depmixS4")


###################################################
### code chunk number 2: depmixS4.Rnw:184-186
###################################################
data("speed")
plot(as.ts(speed[1:168,]), main = "Speed-accuracy trade-off")


###################################################
### code chunk number 3: depmixS4.Rnw:464-469
###################################################
library("depmixS4")
data("speed")
set.seed(1)
mod <- depmix(response = rt ~ 1, data = speed, nstates = 2,
  trstart = runif(4))


###################################################
### code chunk number 4: depmixS4.Rnw:509-510
###################################################
fm <- fit(mod, emc=em.control(rand=FALSE))


###################################################
### code chunk number 5: depmixS4.Rnw:520-521
###################################################
fm 


###################################################
### code chunk number 6: depmixS4.Rnw:537-538
###################################################
summary(fm)


###################################################
### code chunk number 7: depmixS4.Rnw:571-575
###################################################
set.seed(1)
mod <- depmix(rt ~ 1, data = speed, nstates = 2, family = gaussian(),
  transition = ~ scale(Pacc), instart = runif(2)) 
fm <- fit(mod, verbose = FALSE, emc=em.control(rand=FALSE)) 


###################################################
### code chunk number 8: depmixS4.Rnw:582-583
###################################################
summary(fm, which = "transition") 


###################################################
### code chunk number 9: depmixS4.Rnw:601-606
###################################################
set.seed(1)
mod <- depmix(list(rt ~ 1,corr ~ 1), data = speed, nstates = 2, 
  family = list(gaussian(), multinomial("identity")),
  transition = ~ scale(Pacc), instart = runif(2))
fm <- fit(mod, verbose = FALSE, emc=em.control(rand=FALSE))


###################################################
### code chunk number 10: depmixS4.Rnw:611-612
###################################################
summary(fm, which = "response")


###################################################
### code chunk number 11: depmixS4.Rnw:662-663
###################################################
setpars(mod, value = 1:npar(mod))


###################################################
### code chunk number 12: depmixS4.Rnw:669-670
###################################################
setpars(mod, getpars(mod, which = "fixed"))


###################################################
### code chunk number 13: depmixS4.Rnw:675-680
###################################################
trst <- c(0.9, 0.1, 0, 0, 0.1, 0.9, 0, 0)
mod <- depmix(list(rt ~ 1,corr ~ 1), data = speed, transition = ~ Pacc,
  nstates = 2, family = list(gaussian(), multinomial("identity")),
  trstart = trst, instart = c(0.99, 0.01))
fm1 <- fit(mod,verbose = FALSE, emc=em.control(rand=FALSE))


###################################################
### code chunk number 14: depmixS4.Rnw:686-695
###################################################
pars <- c(unlist(getpars(fm1)))
pars[6] <- pars[10] <- 11
pars[1] <- 0
pars[2] <- 1
pars[13] <- pars[14] <- 0.5
fm1 <- setpars(mod, pars)
conpat <- c(0, 0, rep(c(0, 1), 4), 1, 1, 0, 0, 1, 1, 1, 1)
conpat[6] <- conpat[10] <- 2
fm2 <- fit(fm1, equal = conpat)


###################################################
### code chunk number 15: depmixS4.Rnw:768-777
###################################################
data("balance")
set.seed(1)
mod <- mix(list(d1 ~ 1, d2 ~ 1, d3 ~ 1, d4 ~ 1), data = balance,
  nstates = 3, family = list(multinomial("identity"),
  multinomial("identity"), multinomial("identity"),
  multinomial("identity")), respstart = runif(24), prior = ~ age,
  initdata = balance)
fm <- fit(mod, verbose = FALSE, emc=em.control(rand=FALSE)) 
fm


###################################################
### code chunk number 16: depmixS4.Rnw:795-796
###################################################
summary(fm, which = "prior")


###################################################
### code chunk number 17: depmixS4.Rnw:814-827
###################################################
x <- mlogit(base=1)
coeff <- coefficients(fm@prior@parameters)

pr1 <- function(y) {sapply(y, function(z) {x$linkinv(c(t(coeff)%*%c(1,z)), base=1)[1]})}
pr2 <- function(y) {sapply(y, function(z) {x$linkinv(c(t(coeff)%*%c(1,z)), base=1)[2]})}
pr3 <- function(y) {sapply(y, function(z) {x$linkinv(c(t(coeff)%*%c(1,z)), base=1)[3]})}

plot(pr1,min(balance$age),max(balance$age),lty=1,ylim=c(0,1),
main="Prior probabilities by age, balance scale data", xlab="age", ylab="Pr")
plot(pr2,min(balance$age),max(balance$age),add=T,lty=2)
plot(pr3,min(balance$age),max(balance$age),add=T,lty=3)

legend("right",legend=c("Class 1 (correct)","Class 2 (incorrect)","Class 3 (guess)"),lty=1:3,inset=c(0.1,0))


###################################################
### code chunk number 18: depmixS4.Rnw:900-901
###################################################
setClass("exgaus", contains="response")


###################################################
### code chunk number 19: depmixS4.Rnw:925-949
###################################################
library("gamlss")
library("gamlss.dist")
setGeneric("exgaus", function(y, pstart = NULL, fixed = NULL, ...) 
  standardGeneric("exgaus"))

setMethod("exgaus", 
  signature(y = "ANY"), 
  function(y, pstart = NULL, fixed = NULL, ...) {
    y <- matrix(y, length(y))
    x <- matrix(1) 
    parameters <- list()
    npar <- 3
    if(is.null(fixed)) fixed <- as.logical(rep(0, npar))
    if(!is.null(pstart)) {
      if(length(pstart) != npar) stop("length of 'pstart' must be ", npar)
      parameters$mu <- pstart[1]
      parameters$sigma <- log(pstart[2])
      parameters$nu <- log(pstart[3])
    }
    mod <- new("exgaus", parameters = parameters, fixed = fixed,
      x = x, y = y, npar = npar)
    mod
  }
)


###################################################
### code chunk number 20: depmixS4.Rnw:952-1004
###################################################
setMethod("dens","exgaus",
    function(object,log=FALSE) {
    	dexGAUS(object@y, mu = predict(object), 
	sigma = exp(object@parameters$sigma), 
	nu = exp(object@parameters$nu), 
	log = log)
    }
)

setMethod("getpars","response",
    function(object,which="pars",...) {
        switch(which,
            "pars" = {
                parameters <- numeric()
                parameters <- unlist(object@parameters)
                pars <- parameters
            },
            "fixed" = {
                pars <- object@fixed
            }
        )
        return(pars)
    }
)

setMethod("setpars","exgaus",
    function(object, values, which="pars", ...) {
        npar <- npar(object)
        if(length(values)!=npar) stop("length of 'values' must be",npar)
        # determine whether parameters or fixed constraints are being set
		nms <- names(object@parameters)
		switch(which,
		  "pars"= {
		      object@parameters$mu <- values[1]
		      object@parameters$sigma <- values[2]
		      object@parameters$nu <- values[3]
		      },
		  "fixed" = {
		      object@fixed <- as.logical(values)
		  }
	    )
        names(object@parameters) <- nms
        return(object)
    }
)

setMethod("predict","exgaus", 
    function(object) {
        ret <- object@parameters$mu
        return(ret)
    }
)


###################################################
### code chunk number 21: depmixS4.Rnw:1009-1024
###################################################
setMethod("fit", "exgaus",
  function(object, w) {
    if(missing(w)) w <- NULL
    y <- object@y
    fit <- gamlss(y ~ 1, weights = w, family = exGAUS(),
      control = gamlss.control(n.cyc = 100, trace = FALSE),
      mu.start = object@parameters$mu,
      sigma.start = exp(object@parameters$sigma),
      nu.start = exp(object@parameters$nu))    
    pars <- c(fit$mu.coefficients, fit$sigma.coefficients, 
      fit$nu.coefficients)
    object <- setpars(object,pars)
    object
  }
)


###################################################
### code chunk number 22: depmixS4.Rnw:1038-1048
###################################################
rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- exgaus(speed$rt, pstart = c(5, 0.1, 0.1))
rModels[[1]][[2]] <- GLMresponse(formula = corr ~ 1, data = speed,
  family = multinomial(), pstart = c(0.5, 0.5))

rModels[[2]] <- list()
rModels[[2]][[1]] <- exgaus(speed$rt, pstart = c(6, 0.1, 0.1))
rModels[[2]][[2]] <- GLMresponse(formula = corr ~ 1, data = speed, 
  family = multinomial(), pstart = c(0.1, 0.9))


###################################################
### code chunk number 23: depmixS4.Rnw:1054-1062
###################################################
trstart <- c(0.9, 0.1, 0.1, 0.9)
transition <- list()
transition[[1]] <- transInit(~ Pacc, nst = 2, data = speed, 
  pstart = c(0.9, 0.1, 0, 0))
transition[[2]] <- transInit(~ Pacc, nst = 2, data = speed, 
  pstart = c(0.1, 0.9, 0, 0))
inMod <- transInit(~ 1, ns = 2, pstart = c(0.1, 0.9),
  data = data.frame(1))


###################################################
### code chunk number 24: depmixS4.Rnw:1067-1070
###################################################
mod <- makeDepmix(response = rModels, transition = transition,
  prior = inMod, homogeneous = FALSE)
fm <- fit(mod, verbose = FALSE, emc=em.control(rand=FALSE))


