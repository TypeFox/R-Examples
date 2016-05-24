#######################################################################
## Copyright 2008, Novartis Pharma AG
##
## This program is Open Source Software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see http://www.gnu.org/licenses/.

# functions to design trials for MCPMod
library(mvtnorm)
library(lattice)

guesst <-
  function(d, p, model = c("emax", "exponential", "logistic", "quadratic",
                   "betaMod", "sigEmax"), less = TRUE,  local = FALSE, dMax, Maxd, scal)
{
  ## get guesstimates for standardized model parameters(s) 
  ## d = dose corresponding to prior expected percentage "p" of maximum effect
  ## p = percentage corresponding to "d"
  ## model = name of dose response model
  ## less = logical only needed in case of quadratic model.
  ##        determines if d is smaller (less=T) or larger (less=F)
  ##        than d_opt
  ## local = logical indicating whether local (ie exact) or asymptotic versions
  ##         of guesstimate should be derived for emax, logistic and
  ##         sigEmax model (for other models method is exact)
  ##         Branson et al. derive the asymptotic
  ##         expressions (ie. assuming that f^0(Maxd)=1)
  ##         in case local version is used, convergence problems
  ##         with underlying nonlinear optimization can occur
  ## dMax = dose at which max effect occurs (needed only for betaMod)
  ## Maxd = Maximum dose administered
  ## scal = dose scale (needed only for betaMod)
  model <- match.arg(model)
  if(any(p <= 0) | any(p > 1)) stop("must have 0 < p <= 1")
  switch(model,
         emax = {
           if(!local){
             c(ed50 = mean(d * (1 - p)/p))
           } else {
             if (any(p <= d/Maxd))
               stop("must have p > d/Maxd, for local version")
             val <- (d/p-d)/(1-d/(Maxd*p))
             c(ed50=mean(val))
            }
         },
         exponential = {
           if(any(p >= d/Maxd)) stop("must have p < d/Maxd")
           init <- d/log(1 + p)
           fooexp <- function(delta, d, p, Maxd){
             sum((exponential(d, 0, 1, delta)/
                 exponential(Maxd, 0, 1, delta) - p)^2)           
           }
           val <- optimize(fooexp, c(0, 2*Maxd), d=d, p=p, Maxd=Maxd)$minimum
           c(delta = mean(val))
         },
         logistic = {                   
           if(length(d) == 1) {
             stop("logistic model needs at least two pairs (d,p)")
           }
           logit <- function(p) log(p/(1-p)) 
           if(length(d) == 2) {
             ed50  <- diff(rev(d)*logit(p))/diff(logit(p))
             delta <- diff(d)/diff(logit(p))
             res <- c(ed50 = ed50, delta = delta)
           } else {
             m <- lm(logit(p)~d)
             par <- coef(m)
             names(par) <- NULL
             res <- c(ed50 = -par[1]/par[2], delta = 1/par[2])
           } 
           if(local){
             foolog <- function(par, d, p, Maxd) {
               e0 <- logistic(0,0,1,par[1],par[2])
               sum(((logistic(d,0,1,par[1],par[2]) - e0)/
                    (logistic(Maxd,0,1,par[1],par[2])-e0)-p)^2)
             }
             res <- try(optim(par=res, fn=foolog, d=d, p=p, Maxd=Maxd))
             if(res$convergence > 0)
               stop("cannot find guesstimates for specified values")
             else res <- res$par
           }
           if(res[1] < 0)
             message("specified values lead to negative ed50, which should be positive")
           res
         },
         quadratic = {
           aux <- sqrt(1 - p)
           if (less) c(delta = mean(-(1 - aux)/(2 * d)))
           else  c(delta = mean(-(1 + aux)/(2 * d)))
         },
         betaMod = {
           if(scal <= dMax)
             stop("scal needs to be larger than dMax to calculate guesstimate")
           k <- dMax/(scal-dMax)
           val <- d^k*(scal-d)/(dMax^k*(scal-dMax))
           beta <- log(p)/(log(val))
           c(delta1 = mean(k*beta), delta2 = mean(beta))
         },
         sigEmax = {
           if(length(d) == 1) {
             stop("sigmoid Emax model needs at least two pairs (d,p)")
           }
           if(length(d) == 2){
             num <- log((p[1]*(1-p[2]))/(p[2]*(1-p[1])))
             h <- num/log(d[1]/d[2])
             ed50 <- ((1-p[1])/p[1])^(1/h)*d[1]
             res <- c(ed50=ed50, h=h)
           } else {
             y <- log((1-p)/p)
             x <- log(d)
             par <- coef(lm(y~x))
             names(par) <- NULL
             res <- c(ed50 = exp(par[1]/-par[2]), delta = -par[2])
           } 
           if(local) {
             fooSE <- function(par, d, p, Maxd) {
               sum((sigEmax(d,0,1,par[1],par[2])/
                    sigEmax(Maxd,0,1,par[1],par[2])-p)^2)
             }
             res <- try(optim(par=res, fn=fooSE, d=d, p=p, Maxd=Maxd))
             if(res$convergence > 0)
               stop("cannot find guesstimates for specified values")
             else res <- res$par
           }
           if(res[1] < 0)
             message("specified values lead to negative ed50, which should be positive")
           res
         })
}
### example call
### guesst(0.5, 0.8, "emax")
### guesst(0.5, 0.8, "emax", Maxd=1, local=T)
### guesst(0.5, 0.4, "exponential", Maxd = 1)
### guesst(0.7, 1, "quadratic") # maximum response at dose 0.7
### guesst(c(0.4, 0.7), c(0.5, 0.95), "logistic")
### guesst(c(0.4, 0.7), c(0.5, 0.95), "logistic", local = T, Maxd = 1)
### guesst(c(0.4, 0.7), c(0.5, 0.95), "sigEmax")
### guesst(c(0.4, 0.7), c(0.5, 0.95), "sigEmax", local = T, Maxd = 1)
### guesst(0.4,0.5, "betaMod", dMax = 0.8, scal=1.3)

linear <-
  function(dose, e0, delta) e0 + delta * dose

linlog <-
  function(dose, e0, delta, off = 1) linear(log(dose + off), e0, delta)

emax <-  function(dose, e0, eMax, ed50){
  sigEmax(dose, e0, eMax, ed50, 1)
}

quadratic <-
  function(dose, e0, b1, b2) e0 + b1 * dose + b2 * dose^2

exponential <- function(dose, e0, e1, delta){
  e0 + e1*(exp(dose/delta) - 1)
}

logistic <-
  function(dose, e0, eMax, ed50, delta)
{ 
  e0 + eMax/(1 + exp((ed50 - dose)/delta))
}

betaMod <-
  function(dose, e0, eMax, delta1, delta2, scal)
{
  maxDens <- (delta1^delta1)*(delta2^delta2)/
    ((delta1 + delta2)^(delta1+delta2))
  dose <- dose/scal
  e0 + eMax/maxDens * (dose^delta1) * (1 - dose)^delta2
}

sigEmax <- 
  function(dose, e0, eMax, ed50, h)
{
  e0 + eMax*dose^h/(ed50^h + dose^h)
}

## means corresponding to different models and doses
modelMeans <-
  ## generate mean vectors for models and initial values defined by user
  function(models, doses, std = TRUE, off = 0.1*max(doses), scal = 1.2*max(doses))
{
  ## models - list with names given by the model functions to be used
  ##          and components as eitheir vectors or matrices with the
  ##          parameter values ; assumption is that first argument to
  ##          model function is always the doses 
  ## doses  - doses to be used in design
  ## std    - logical variable indicating whether to assume
  ##          standardized versions of built-in models
  
  namMod <- names(models)
  if(length(namMod) != length(unique(namMod))){
    stop("only one list entry allowed for each model class in 'models' argument")
  }
  ## built-in models with methods for standardized versions
  biModels <- c("emax", "linlog", "linear", "quadratic",
                 "exponential", "logistic", "betaMod", "sigEmax")
  nModels <- length(models)             # number of model elements
  contMap <- vector("list", nModels)
  names(contMap) <- namMod
  val <- list()
  k <- 1
  nams <- character()
  for(nm in namMod) {
    pars <- models[[nm]]
    if (!is.null(pars) && !is.numeric(pars)) {
      stop("elements of \"models\" must be NULL or numeric")
    }
    if (is.element(nm, biModels) && std) {
      if (!is.element(nm, c("betaMod","logistic","sigEmax"))) {
        ## all others have single std parameter
        if (is.matrix(pars)) {
          stop("For standardized ", nm,
                     " models element cannot be matrix")
        }
        nmod <- length(pars)
        if (nmod > 1) {                 # multiple models
          ind <- 1:nmod
          nams <- c(nams, paste(nm, ind, sep = ""))
          contMap[[nm]] <- (k - 1) + ind
          for(j in 1:nmod) {
            if (nm == "linlog") pars1 <- c(0, 1, off)
            else pars1 <- c(0, 1, pars[j])
            val[[k]] <- do.call(nm, c(list(doses), as.list(pars1)))
            k <- k + 1
          }
        } else {                        # single model
          nams <- c(nams, nm)
          contMap[[nm]] <- k
          if (nm == "quadratic") names(pars) <- NULL
          if (nm == "linlog") pars <- c(0, 1, off)
          else pars <- c(0, 1, pars)
          val[[k]] <- do.call(nm, c(list(doses), as.list(pars)))
          k <- k + 1
        }
      } else {                          # logistic, betaMod
        if (is.matrix(pars)) {          # multiple models
          if (ncol(pars) != 2) {
            stop("standardized ", nm," model must have two parameters")
          }
          nmod <- nrow(pars)            # number of models
          ind <- 1:nmod
          nams <- c(nams, paste(nm, ind, sep = ""))
          contMap[[nm]] <- (k - 1) + ind
          for(j in 1:nmod) {
            if(nm == "betaMod")     pars1 <- c(0, 1, pars[j, ], scal)
            else pars1 <- c(0, 1, pars[j, ])
            val[[k]] <- do.call(nm, c(list(doses), as.list(pars1)))
            k <- k + 1
          }
        } else {                        # single model
          if (length(pars) != 2) {
            stop("standardized ", nm," model must have two parameters")
          }
          nams <- c(nams, nm)
          contMap[[nm]] <- k
          if(nm == "betaMod") pars <- c(pars, scal)
          val[[k]] <-
            do.call(nm, c(list(doses), as.list(c(0, 1, pars))))
          k <- k + 1
        }
      }
    } else {                        # user-defined or non-standardized
      if (is.matrix(pars)) {            # multiple models
        nmod <- nrow(pars)              # number of models
        if(nm == "linlog")  pars <- cbind(pars, off)
        if(nm == "betaMod") pars <- cbind(pars, c(scal))
        ind <- 1:nmod
        nams <- c(nams, paste(nm, ind, sep = ""))
        contMap[[nm]] <- (k - 1) + ind
        for(j in 1:nmod) {
          val[[k]] <- do.call(nm, c(list(doses), as.list(pars[j,])))
          k <- k + 1
        }
      } else {                      # single model
        if(nm == "linlog")  pars <- c(pars, off)
        if(nm == "betaMod") pars <- c(pars, scal)                     
        nams <- c(nams, nm)
        contMap[[nm]] <- k
        val[[k]] <- do.call(nm, c(list(doses), as.list(pars)))
        k <- k + 1
      }       
    }
  }
  muMat <- do.call("cbind", val)
  dimnames(muMat) <- list(doses, nams)
  attr(muMat, "contMap") <- contMap
  muMat
}

## Model contrasts
getStdCont <-
  function(cont)
{
  ## center and normalize vector
  cont <- cont - mean(cont)
  cont/sqrt(sum(cont^2))
}

### simplified version
getOptCont <-
  function(mu, n = 1)
{
  mn <- sum(mu * n)/sum(n)
  getStdCont(n * (mu - mn))
}

modContr <-
  function(means, n = 1)
{
  ## obtains model contrasts corresponding to model means and
  ## sample sizes n
  apply(means, 2, getOptCont, n = n)
}

# control parameters for functions in package mvtnorm
mvtnorm.control <-
  function(maxpts = 30000, abseps = 0.001,
           interval = c(0,10), releps = 0)
{
  out <- list(interval = interval, maxpts = maxpts, abseps = abseps,
       releps = releps)
  class(out) <- "GenzBretz"
  out
}
  
critVal <-
  function(cMat, n, alpha = 0.025, control = mvtnorm.control(),
           twoSide = FALSE, corMat = NULL)
{
  ## function to calculate critical value
  ##
  ## cMat - model contrast matrix nDose x nMod
  ## n - scalar or vector with sample size(s)
  ## alpha - significance level
  ## control - see mvtnorm.control function
  ## twoSide - logical variable indicating if one or two-sided test
  ## corMat - correlation matrix

  nD <- nrow(cMat)
  nMod <- ncol(cMat)
  if (length(n) == 1) {
    nDF <- nD * (n - 1)
  } else {
    if (length(n) != nD) {
      stop("sample size vector must have length equal to number of doses")
    }
    nDF <- sum(n) - nD
  }
  if(is.null(corMat)){
      corMat <- t(cMat)%*%(cMat/n)
      den  <- sqrt(crossprod(t(colSums(cMat^2/n))))
      corMat <- corMat / den
  }
  if(twoSide) tail <- "both.tails"
  else tail <- "lower.tail"
  if (!missing(control)) {
    if(!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  } else {
    ctrl <- control
  }
  qmvtCall <- c(list(1-alpha, tail = tail, df = nDF, corr = corMat, 
                     algorithm=ctrl, interval=ctrl$interval))
  do.call("qmvt", qmvtCall)$quantile
}

planMM <-
  function(models, doses, n, off = 0.1*max(doses), scal = 1.2*max(doses), std = TRUE,
           alpha = 0.025, twoSide = FALSE, control = mvtnorm.control(), 
           cV = TRUE, muMat = NULL)
{
  ## returns the critical value, the contrast and correlation matrices
  ## models - list with names given by the model functions to be used
  ##          and components as either vectors or matrices with the
  ##          parameter values ; assumption is that first argument to
  ##          model function is always the doses 
  ## doses  - doses to be used in design
  ## n - scalar or vector with sample size(s)
  ## off, scal - additional parameters for the linear in log and beta model
  ## std    - logical variable indicating whether to assume
  ##          standardized versions of built-in models
  ## alpha - significance level
  ## twoSide - logical variable indicating if one or two-sided test
  ## control - control arguments for mvtnorm function(s)
  ## cV - logical indicating whether critical value should be calc.
  ## muMat - an optional matrix with means as columns, dimnames should
  ##         be given (dose levels and names of contrasts)

  if (missing(models)) {
    ## may pass necessary information via muMat
    if (is.null(muMat)) {
      stop("either models or muMat need to be specified")
    }
    mu <- muMat
  } else {
    mu <- modelMeans(models, doses, std, off, scal)
  }
  contMat <- modContr(mu, n)
  corMat <- t(contMat)%*%(contMat/n)
  den  <- sqrt(crossprod(t(colSums(contMat^2/n))))
  corMat <- corMat / den
  if(cV){
    critV <- critVal(contMat, n, alpha, control = control, twoSide = twoSide)
  }
  else critV <- NULL
  res <- list(contMat = contMat, critVal = critV, muMat = mu,
              corMat = (corMat+t(corMat))/2 )
  attr(res, "n") <- n
  attr(res, "alpha") <- alpha
  if (twoSide) {
    attr(res, "side") <- "two-sided"
  } else {
    attr(res, "side") <- "one-sided"
  }
  oldClass(res) <- "planMM"
  res
}

print.planMM <-
  function(x, digits = 3, ...)
{
  cat("MCPMod planMM\n")
  cat("\n","Optimal Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  if(!is.null(x$critVal)){
    names(x$critVal) <- NULL
    cat("\n","Critical Value ", sep="")
    cat(paste("(alpha = ", attr(x, "alpha"),", ", attr(x, "side"), "): ",
              round(x$critVal, digits), sep=""),"\n")
  }
  cat("\n","Contrast Correlation Matrix:","\n", sep="")
  print(round(x$corMat, digits))  
  cat("\n")
}

### example call (example from JBS paper)
### doses <- c(0, 10, 25, 50, 100, 150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential = c(85),
###                betaMod = matrix(c(0.33, 2.31, 1.39, 1.39), byrow=T,nrow=2))
### planMM(models, doses, 50, scal = 200, alpha = 0.05) # compare with JBS 16, 645-646
### example with muMat         
### dvec <- c(0, 10, 50, 100)
### mu1 <- c(1, 2, 2, 2)
### mu2 <- c(1, 1, 2, 2)
### mu3 <- c(1, 1, 1, 2)
### mMat <- cbind(mu1, mu2, mu3)
### dimnames(mMat)[[1]] <- dvec
### planMM(muMat = mMat, doses = dvec, n = 30)

plot.planMM <- function (x, superpose = TRUE, xlab = "Dose",
    ylab = NULL, resp = c("contrasts", "means"), ...) 
{
    resp <- match.arg(resp)
    if (is.null(ylab)) {
        if (resp == "contrasts") {
            ylab <- "Contrast coefficients"
        }
        else {
            ylab <- "Normalized model means"
        }
    }
    cM <- x$contMat
    if (resp == "means") {
        cM <- t(t(x$muMat)/apply(x$muMat, 2, max))
    }
    nD <- nrow(cM)
    nM <- ncol(cM)
    cMtr <- data.frame(resp = as.vector(cM), dose = rep(as.numeric(dimnames(cM)[[1]]), 
        nM), model = factor(rep(dimnames(cM)[[2]], each = nD), 
        levels = dimnames(cM)[[2]]))
    if (superpose) {
        spL <- trellis.par.get("superpose.line")
        spL$lty <- rep(spL$lty, nM%/%length(spL$lty) + 1)[1:nM]
        spL$lwd <- rep(spL$lwd, nM%/%length(spL$lwd) + 1)[1:nM]
        spL$col <- rep(spL$col, nM%/%length(spL$col) + 1)[1:nM]
        xyplot(resp ~ dose, data = cMtr, subscripts = TRUE, groups = cMtr$model, 
            panel = panel.superpose, type = "o", xlab = xlab, 
            ylab = ylab, key = list(lines = spL, transparent = TRUE, 
                text = list(levels(cMtr$model), cex = 0.9), columns = ifelse(nM < 
                  5, nM, min(4,ceiling(nM/min(ceiling(nM/4),3))))), ...)
    }
    else {
        xyplot(resp ~ dose | model, data = cMtr, type = "o", 
            xlab = xlab, ylab = ylab, strip = function(...) strip.default(..., 
                style = 1), ...)
    }
}

### doses <- c(0, 10, 25, 50, 100, 150)                                        
### models <- list(linear = NULL, emax = c(25),                                
###                logistic = c(50, 10.88111), exponential = c(85),            
###                betaMod = matrix(c(0.33, 2.31, 1.39, 1.39), byrow=T,nrow=2))
### pM <- planMM(models, doses, 50, scal = 200)
### plot(pM)
### plot(pM, superpose=F, xlab="Different axis name")
### plot(pM, resp = "means")
### #example with muMat
### dvec <- c(0, 10, 50, 100)
### mu1 <- c(1, 2, 2, 2)
### mu2 <- c(1, 1, 2, 2)
### mu3 <- c(1, 1, 1, 2)
### mMat <- cbind(mu1, mu2, mu3)
### dimnames(mMat)[[1]] <- dvec
### pM <- planMM(muMat = mMat, doses = dvec, n = 30)
### plot(pM)                                         
### plot(pM, superpose=F, xlab="Different axis name")
### plot(pM, resp = "means")                         

getPars <-
  function(model, doses, initEstim, base, maxEff, off = 0.1*max(doses), scal = 1.2*max(doses))
{
  ## calculates the location and scale parameters corresponding to
  ## given base, maxEff, and guesstimates
  ## model - a character string with the model name. Built-in models
  ##         have their full parameterization derived internally. For
  ##         user-defined models, it is assumed the existence of a
  ##         function named as "Par" appended to end of model (e.g.,
  ##         for model = "cubic", it is assumed that there is function
  ##         cubicPar that calculates the necessary parameters; this
  ##         function is assumed to have arguments doses, initEstim,
  ##         base, and maxEff, in that order. (see below for an example)
  ## doses - doses to be used in design
  ## base, maxEff - baseline and maximum effect (over placebo)
  ## initEstim - vector of guesstimates
  ## off, scal - additional parameters for the linear in log and beta model
  switch(model,
         linear = {
           e1 <- maxEff/max(doses)
           Par <- c(base, e1)
         },
         linlog = {
           e1 <- maxEff/(log(max(doses) + off) - log(off))
           Par <- c(base-e1*log(off), e1)
         },
         quadratic = {
           dMax <- 1/(-2*initEstim)
           b1 <- maxEff/(dMax + initEstim*dMax^2)
           b2 <- initEstim * b1
           Par <- c(base, b1, b2)
         },
         emax = {
           emax.p <- maxEff * (initEstim + max(doses))/max(doses)
           Par <- c(base, emax.p, initEstim)
         },
         exponential = {
           e1 <- maxEff/(exp(max(doses)/initEstim) - 1)
           e0 <- base
           Par <- c(e0, e1, initEstim)
         },
         logistic = {
           emax.p <- maxEff/
             (logistic(max(doses),0,1, initEstim[1], initEstim[2]) -
              logistic(0, 0, 1, initEstim[1], initEstim[2]))
           e0 <- base-emax.p*logistic(0,0,1,initEstim[1], initEstim[2])
           Par <- c(e0, emax.p, initEstim[1], initEstim[2])
         },
         betaMod = {
           Par <- c(base, maxEff, initEstim)
         },
         sigEmax = {
           ed50 <- initEstim[1]
           h <- initEstim[2]
           dmax <- max(doses)
           eMax <- maxEff*(ed50^h+dmax^h)/dmax^h
           Par <- c(e0 = base, eMax = eMax, ed50 = ed50, h = h)
         },
         {
           Par <- do.call(paste(model,"Par", sep=""),
                          list(doses, initEstim, base, maxEff))
         }
         )
  Par
}

### example call
### doses <- c(0, 10, 25, 50, 100, 150)
### getPars("emax", doses, 25, 0, 0.4)
### getPars("logistic", doses, c(50, 10.88111), 0, 0.4) # compare JBS 16, p.650
### getPars("betaMod", doses, initEstim = c(0.33, 2.31), base = 0,
###         maxEff = 0.4)
### #example for user model
### sigEmax2 <- function(dose, e0, eMax, ed50, h){
###   e0 + eMax * ( dose^h / (ed50^h + dose^h) )
### }
### # function to return location and scale parameters
### sigEmax2Par <- function(dose, initEstim, base, maxEff){
###   # function to get linear parameters
###   # ed50 parameter assumed to be first in initEstim
###   ed50 <- initEstim[1]
###   h <- initEstim[2]
###   dmax <- max(dose)
###   emax <- maxEff*(ed50^h+dmax^h)/dmax^h
###   c(base, emax, initEstim)
### }
### getPars("sigEmax2", doses, initEstim = c(50,2), base = 0, maxEff = 1)


fullMod <-
  function(models, doses, base, maxEff, off = 0.1*max(doses), scal = 1.2*max(doses))
{
  ## calculates parameters for all models in the candidate set
  ## returns a list (similar to the models list) but with all model parameters.
  ## models - list with names given by the model functions to be used
  ##          and components as either vectors or matrices with the
  ##          parameter values
  ## doses  - doses to be used in design
  ## base, maxEff - baseline, maximum effect (over placebo)
  ## off, scal - additional parameters for the linear in log and beta model
  
  namMod <- names(models)
  if(length(namMod) != length(unique(namMod))){
    stop("only one list entry allowed for each model class in 'models' argument")
  }  
  complModels <- list()
  i <- 0
  for(nm in namMod){
    pars <- models[[nm]]
    if(is.null(pars)){
      Pars <- getPars(nm, doses, NULL, base, maxEff, off); i <- i+1
    }
    else{ 
      if(!is.element(nm,c("logistic", "betaMod", "sigEmax"))){
        nmod <- length(pars)
        if(nmod > 1){
          Pars <- matrix(ncol=3, nrow=nmod)
          for(j in 1:length(pars)){
            Pars[j,] <-  getPars(nm, doses, as.vector(pars[j]), base, maxEff)
          }
          i <- i+1
        }
        else {
          Pars <-  getPars(nm, doses, as.vector(pars), base, maxEff)
          i <- i+1
        }
      }
      else { 
        if(is.matrix(pars)){
          Pars <- matrix(ncol=4, nrow=nrow(pars))
          for(j in 1:nrow(pars)){
            Pars[j,] <-  getPars(nm, doses, as.vector(pars[j,]), base, maxEff)
          }
          i <- i+1
        }
        else {
          Pars <-  getPars(nm, doses, as.vector(pars), base, maxEff); i <- i+1
        }
      }
    }
    complModels[[i]] <- Pars
  }
  names(complModels) <- names(models)
  ## store information for use later
  attr(complModels, "doses") <- doses
  attr(complModels, "base") <- base
  attr(complModels, "maxEff") <- maxEff
  attr(complModels, "off") <- off
  attr(complModels, "scal") <- scal
  oldClass(complModels) <- "fullMod"
  complModels
}

### example call
### doses <- c(0, 10, 25, 50, 100, 150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential = c(85),
###                betaMod = matrix(c(0.33, 2.31, 1.39, 1.39), byrow=T,nrow=2))
### fullMod(models, doses, base = 0, maxEff = 0.4, scal = 200)

plotModels <- function (models, doses, base, maxEff, nPoints = 200, off = 0.1*max(doses), 
    scal = 1.2 * max(doses), superpose = FALSE, ylab = "Model means", 
    xlab = "Dose", ...) 
{
    if (inherits(models, "fullMod")) {
        doses <- attr(models, "doses")
        base <- attr(models, "base")
        maxEff <- attr(models, "maxEff")
        off <- attr(models, "off")
        scal <- attr(models, "scal")
        Models <- models
    }
    else {
        Models <- fullMod(models, doses, base, maxEff, off, scal)
    }
    doseSeq <- sort(union(seq(min(doses), max(doses), length = nPoints), 
        doses))
    resp <- modelMeans(Models, doseSeq, std = FALSE, off, scal)
    nams <- dimnames(resp)[[2]]
    nM <- length(nams)
    if (superpose) {
       respdata <- data.frame(response = c(resp),
                         dose = rep(doseSeq, ncol(resp)),
                         model = factor(rep(nams, each = length(doseSeq)),
                                        levels = nams))
        spL <- trellis.par.get("superpose.line")
        spL$lty <- rep(spL$lty, nM%/%length(spL$lty) + 1)[1:nM]
        spL$lwd <- rep(spL$lwd, nM%/%length(spL$lwd) + 1)[1:nM]
        spL$col <- rep(spL$col, nM%/%length(spL$col) + 1)[1:nM]
        xyplot(response ~ dose, data = respdata, subscripts = TRUE, 
            groups = respdata$model, panel.data = list(base = base, maxEff = maxEff, 
                doses = doses), xlab = xlab, ylab = ylab, panel = function(x, 
                y, subscripts, groups, ..., panel.data) {
                panel.abline(h = c(panel.data$base, panel.data$base + 
                  panel.data$maxEff), lty = 2)
                panel.superpose(x, y, subscripts, groups, type = "l", 
                  ...)
                ind <- !is.na(match(x, panel.data$doses))
                panel.superpose(x[ind], y[ind], subscripts[ind], 
                  groups, ...)
            }, key = list(lines = spL, transparent = TRUE, text = list(nams, 
                cex = 0.9), columns = ifelse(nM < 
                  5, nM, min(4,ceiling(nM/min(ceiling(nM/4),3))))), ...)
    }
    else {
        respdata <- data.frame(response = c(resp), dose = rep(doseSeq, 
              ncol(resp)), model = factor(rep(nams, each = length(doseSeq))))
        xyplot(response ~ dose | model, data = respdata, panel.data = list(base = base, 
            maxEff = maxEff, doses = doses), xlab = xlab, ylab = ylab, 
            panel = function(x, y, ..., panel.data) {
                panel.abline(h = c(panel.data$base, panel.data$base + 
                  panel.data$maxEff), lty = 2)
                panel.xyplot(x, y, type = "l", ...)
                ind <- match(panel.data$doses, x)
                panel.xyplot(x[ind], y[ind], ...)
            }, strip = function(...) strip.default(..., style = 1), 
            as.table = TRUE,...)
    }
}

plot.fullMod <-
  function(x, ...)
{
  plotModels(x, ...)
}

### example call
### doses <- c(0,10,25,50,100,150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential = c(85),
###                betaMod = matrix(c(0.33, 2.31, 1.39, 1.39), byrow=T,nrow=2))
### plotModels(models, doses, base = 0, maxEff = 0.4, scal = 200)
### plotModels(models, doses, base = 0, maxEff = 0.4, scal = 200, superpose = T)

powCalc <-
  function(cMat, n, alpha = 0.025, delta = NULL, mu = NULL,
           sigma = NULL, cVal = NULL, corMat = NULL,
           twoSide = FALSE, control = mvtnorm.control())
{
  ## cMat - contrast matrix with rows (as number of doses) and columns
  ##        (as number of models) - exactly as in the output of planMM
  ## n - scalar or vector with sample size(s)
  ## alpha - significance level
  ## delta - non-centrality vector (in same order as models)
  ## mu - vector of model means (either delta or mu and sigma need to be given)
  ## sigma - standard deviation of response (to build delta vector)
  ## cVal - optional critical value for the test
  ## twoSide - logical variable indicating if one or two-sided test

  nD <- nrow(cMat)
  nMod <- ncol(cMat)
  if (length(n) == 1) {
    nDF <- nD * (n - 1)
  } else {
    if (length(n) != nD) {
      stop("sample size vector must have length equal to number of doses")
    }
    nDF <- sum(n) - nD
  }
  if (is.null(delta)) {    
    den  <- sqrt(colSums(cMat^2/n))
    delta <- as.vector((t(cMat) %*% mu)/(den*sigma))
  }
  if(is.null(corMat)){
      corMat <- t(cMat)%*%(cMat/n)
      den  <- sqrt(crossprod(t(colSums(cMat^2/n))))
      corMat <- corMat / den
  }
  if (is.null(cVal)) {
    cVal <-
      critVal(cMat, n, alpha, control, twoSide, corMat)
  }          
  if(twoSide){
    lower <- rep(-cVal, nMod)
  }
  else lower <- rep(-Inf, nMod)
  upper <- rep(cVal, nMod)
  if (!missing(control)) {
    if(!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  } else {
    ctrl <- control
  }
  ctrl$interval <- NULL      # not used with pmvt
  pmvtCall <- c(list(lower, upper, df = nDF, corr = corMat,
                     delta = delta, algorithm=ctrl))
  as.vector(1 - do.call("pmvt", pmvtCall))
}
### example call
### doses <- c(0,10,25,50,100,150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential=c(85),
###                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=T, nrow=2))
### plMM <- planMM(models, doses, 50, scal = 200, alpha = 0.05)
### compMod <- fullMod(models, doses, base = 0, maxEff = 0.4, scal = 200)
### muMat <- modelMeans(compMod, doses, F, scal = 200)
### powCalc(plMM$contMat, 50, mu = muMat[,3], sigma = 1, cVal = plMM$critVal, control = list(maxpts=100000)) # Power for logistic model
### powCalc(plMM$contMat, 50, mu = muMat[,1], sigma = 1, cVal = plMM$critVal) # Power for linear model
### # compare with JBS 16, 650

powerMM <-
  function(models, doses, base, maxEff, sigma, lower, upper,
           step, sumFct = c("min", "mean", "max"), off = 0.1*max(doses),
           scal = 1.2*max(doses), alpha = 0.025, twoSide = FALSE,
           control = mvtnorm.control(), muMat = NULL, alRatio = NULL,
           typeN = c("arm", "total"), ...)
{
  ## generates a matrix with power values for different sample sizes
  ## (balanced case only) and models
  ## models - either a list with model parameters as in planMM
  ##          or an object of class fullMod
  ## doses, maxEff, base, off, scal - as in getPars function
  ## sigma - expected standard deviation
  ## lower, upper - maximum and minimum group sample size for which the
  ##                power is plotted
  ## step - in which steps should the power be calculated
  ##        (it is calculated at seq(lower,upper,by=step))
  ## sumFct - a character vector giving the names of the summary
  ##          functions to be called (each function should return a numeric
  ##          of length 1)
  ## off,scal - additional parameters for linlog/beta model
  ## alpha - level of significance
  ## twoSide - logical variable indicating if one or two-sided test
  ## control - control list for mvtnorm fcts., see mvtnorm.control()
  ## muMat - an optional matrix with means as columns
  ## alRatio - optional vector with sample size allocation ratios
  ##           must all be > 0; lower and upper are interpreted as
  ##           limit values for dose(s) with smallest alloc. ratios
  ## ... - possible additional arguments for sumFct 

  if (!missing(models)) {
    ## if models is of class fullMod, will ignore other arguments
    if (inherits(models, "fullMod")) {
      doses <- attr(models, "doses")
      base <- attr(models, "base")
      maxEff <- attr(models, "maxEff")
      off <- attr(models, "off")
      scal <- attr(models, "scal")
      Models <- models
    } else {
      Models <- fullMod(models, doses, base, maxEff, off, scal)
    }
    muMat <- modelMeans(Models, doses, FALSE, off, scal)
  } else {
    if (is.null(muMat)) {
      stop("muMat must be specified, when models is not")
    }
  }
  nD <- length(doses)
  # allocation ratios
  Ntype <- match.arg(typeN)
  if(!is.null(alRatio)){
    if(length(alRatio)!= nrow(muMat)){
      stop("alRatio needs to be of same length as number of dose levels")
    }
    if(any(alRatio <= 0)){
       stop("all entries of alRatio need to be positive")
    }
    if (Ntype == "arm") {
       alRatio <- alRatio/min(alRatio)
    } else {
       alRatio <- alRatio/sum(alRatio)
    }
  } else {
    alRatio <- 1
  }
 
  contMat <- modContr(muMat, upper*alRatio)
  nmod <- ncol(contMat)
  pow <- numeric(nmod)
  j <- 1
  nSeq <- seq(lower, upper, by = step)
  powMat <- matrix(nrow=length(nSeq), ncol=nmod)  
  for(n in nSeq) {
    N <- round(n * alRatio)             # simple rounding
    cVal <- critVal(contMat, N, alpha, control, twoSide)
    for(i in 1:nmod) {
      den  <- sqrt(colSums(contMat^2/N))
      delta <- as.vector((t(contMat) %*% muMat[,i])/(den*sigma))
      pow[i] <- powCalc(contMat, N, delta = delta, cVal = cVal,
                        twoSide = twoSide, control = control)
    }
    powMat[j,] <- pow
    j <- j + 1
  }
  nams <- dimnames(muMat)[[2]]
  
  if(!is.character(sumFct)){
    if(length(sumFct) > 1){
      stop("sumFct needs to be a character vector")
    } else {
      sumFct <- deparse(substitute(sumFct)) # for back-compatibility
    }
  }
  nSfunc <- length(sumFct)
  for(i in 1:nSfunc){
    powMat <- cbind(powMat, apply(powMat[,1:nmod, drop = FALSE], 1, sumFct[i] ,...))
  } 

  dimnames(powMat) <- list(nSeq, c(nams, sumFct))
  attr(powMat, "alRatio") <- alRatio
  attr(powMat, "sumFct") <- sumFct
  oldClass(powMat) <- "powerMM"
  powMat
}

### doses <- c(0,10,25,50,100,150)
### models <- list(linear = NULL, emax = 25,
###                logistic = c(50, 10.88111), exponential= 85,
###                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=T, nrow=2))
### pM <- powerMM(models, doses, base = 0, maxEff = 0.4, sigma = 1, alpha = 0.05,
###                lower = 10, upper = 100, step = 10, scal = 200)    
### plot(pM, models="none", line.at = 0.8)

plot.powerMM <-
  function(x, superpose = TRUE, line.at = NULL, models = "all",
           summ = NULL, perc = FALSE, xlab = NULL,
           ylab = ifelse(perc, "Power (%)", "Power"), ...)
{
  ## Trellis plot of power matrix produced by powerMM
  ## x - matrix with power values for different sample sizes and models
  ## superpose - logical variable indicating if lines should be superpose
  ## line.at - a value, or a vector of values, between 0 and 1, to be
  ##           drawn as horizontal line on plot; maybe most interesting if it 
  ##           is set to the desired power value(s) (default: not drawn)
  ## models - which of available models should be included in plot
  ##          "all" and "none" are accepted, else names (or numbers) of models
  ## summ  - summaries to be included in plot; "summary" can be used for
  ##         summary function in x - min and max are also accepted (and
  ##         in combination)
  ## perc - logical indicating if power values should be in percentage
  ## xlab, ylab - labels for x- and y-axis

  ## checking for unbalanced sample sizes
  unbN <- (length(unique(attr(x, "alRatio"))) > 1)
  if (is.null(xlab)) {
    if (unbN) {
      if(sum(attr(x, "alRatio"))==1)
        xlab <- "Overall sample size (unbalanced allocation)"
      else
        xlab <- "Sample size for smallest arm (unbalanced)"
    } else {
      xlab <- "Sample size per dose (balanced)"
    }
  }
  if (perc) x <- 100 * x
  nSeq <- as.integer(dimnames(x)[[1]])
  nams <- dimnames(x)[[2]]
  ## separating model data from summary data
  sumFct <- attr(x, "sumFct")      # original summary functions
  ## default is to use same summary functions as in x
  if (is.null(summ)) summ <- sumFct
  if (!is.null(sumFct)) {
    ## summary functions available in data and requested in summ
    indS <- !is.na(match(sumFct, summ))
    if (any(indS)) {
       summData <- x[, sumFct[indS], drop = FALSE]
    } else {
       summData <- NULL
    }
    modData <- x[, is.na(match(nams, sumFct)), drop = FALSE]
  } else {
    summData <- NULL
    modData <- x
  }
  namsMod <- dimnames(modData)[[2]]
  ## additional summary functions possibly requested in summ
  indSX <- is.na(match(summ, sumFct))
  if(any(indSX)) {
     summ <- summ[indSX]
     if((length(summ) == 1) && (summ == "none")) {
        summDataX <- NULL
     } else {
        nSX <- length(summ)
        summDataX <- vector("list", nSX)
        names(summDataX) <- summ
        for(i in 1:nSX) {
           ## should allow additional arguments to summ[i], but too messy already
           summDataX[[i]] <- apply(modData, 1, summ[i])
        }
        summDataX <- data.frame(summDataX)
     }   
  } else {
     summDataX <- NULL
  }
  ## combining the two summary datasets
  if(length(summData) > 0){
    summData <- cbind(summData, summDataX)
  } else {
    summData <- summDataX
  }
  ## checking subset of models
  if (length(models) == 1) {
    if (models != "all") {
      if (models == "none") {
        modData <- NULL
      } else {
        if (is.numeric(models)) {
          if (!is.element(models, 1:ncol(modData))) {
            stop("invalid model selected")
          }
        } else {
          if (!is.element(models, namsMod)) {
            stop("invalid model selected")
          }
        }
        modData <- modData[,models, drop = FALSE] 
      }
    }
  } else {                              # need to be subset of models
    if (is.numeric(models)) {
      if (!all(is.element(models, 1:ncol(modData)))) {
        stop("invalid models selected")
      }
    } else {
      if (!all(is.element(models, namsMod))) {
        stop("invalid model selected")
      }
    }
    modData <- modData[,models, drop = FALSE]
  }
  if (length(modData) == 0) x <- summData
  else if (length(summData) == 0) x <- modData
  else x <- cbind(modData, summData)
  x <- data.frame(x)
  nams <- names(x)
  nC <- ncol(x)
  pMatTr <-
    data.frame(pow = as.vector(unlist(x)),
               n = rep(nSeq, nC),
               type = factor(rep(nams, each = length(nSeq)), levels = nams))
  if (superpose) {
      panelFunc <- if (is.null(line.at)) {
        panel.superpose
      } else {
        function(x, y, subscripts, groups, lineAt, ...) {
          panel.superpose(x, y, subscripts, groups, ...)
          panel.abline(h = lineAt, lty = 2)
        }
    }
    trLn <- trellis.par.get("superpose.line")[c("col", "lwd", "lty")]
    for(i in seq(along = trLn)) {
       if(length(trLn[[i]]) > nC) trLn[[i]] <- trLn[[i]][1:nC]
    }
    xyplot(pow ~ n, pMatTr, groups = pMatTr$type, subscripts = TRUE,
           panel = panelFunc, type = "l", lineAt = line.at,
           xlab = xlab, ylab = ylab,
           key = list(lines = trLn, text = list(lab = nams), transparent = TRUE, 
           columns = ifelse(nC < 5, nC, min(4,ceiling(nC/min(ceiling(nC/4),3))))), ...)
  } else {                              # models in different panels
    if (is.null(line.at)) {
      panelFunc <- panel.xyplot
    }
    else {
      panelFunc <-
        function(x, y, lineAt, ...) {
          panel.xyplot(x, y, ...)
          panel.abline(h = lineAt, lty = 2) ## used 2 for consistency with above
        }
    }
    xyplot(pow ~ n | type, pMatTr, panel = panelFunc,
           type = "l", lineAt = line.at,
           xlab = xlab, ylab = ylab, 
           strip = function(...) strip.default(..., style = 1), ...)
  }
}

### example call
### doses <- c(0,10,25,50,100,150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential=c(85),
###                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=T, nrow=2))
### pM <- powerMM(models, doses, base = 0, maxEff = 0.4, sigma = 1,
###                lower = 10, upper = 100, step = 5, scal = 200)
### plot(pM)         
### # reproduces plot in JBS 16, p.651
### plot(pM, line.at = 0.8, models = "none") 

sampSize <-
  function(models, doses, base, maxEff, sigma, upperN, lowerN = floor(upperN/2),
           power = 0.8, alRatio = NULL, sumFct = mean, off = 0.1*max(doses), scal = 1.2*max(doses),
           alpha = 0.025, twoSide = FALSE, tol = 1e-3, verbose = FALSE, 
           control = mvtnorm.control(), muMat = NULL, typeN = c("arm", "total"), ...)
{
  ## Calculates  smallest sample size necessary to achieve given power
  ## using a bisection search
  ## models - either a list with model parameters as in planMM
  ##          or an object of class fullMod
  ## doses, maxEff, base, off, scal - as in getPars function
  ## sigma - expected standard deviation
  ## lowerN, upperN - lower and upper limit for the necessary sample size
  ## power - desired power level (between 0 and 1)
  ## alRatio - allocation ratios ()
  ## sumFct - summary function for the individual power values
  ## alpha - level of significance
  ## twoSide - logical variable indicating if one or two-sided test
  ## tol - tolerance for power calculation
  ## verbose - logical controlling amount of output in iterations
  ## control - see mvtnorm.control function
  ## typeN - type of sample size specified in lowerN and upperN, in the case of
  ##         unequal allocation: "arm" refers to smallest trt arm, "total" to
  ##         overall sample size
  ## ... - possible additional arguments for sumFct 
  
  if (missing(models)) {
    ## may pass necessary information via muMat
    if (is.null(muMat)) {
      stop("either models or muMat need to be specified")
    }
  } else {
    ## if models is of class fullMod, will ignore other arguments
    if (inherits(models, "fullMod")) {
       doses <- attr(models, "doses")
       base <- attr(models, "base")
       maxEff <- attr(models, "maxEff")
       off <- attr(models, "off")
       scal <- attr(models, "scal")
       Models <- models
    } else {
       Models <- fullMod(models, doses, base, maxEff, off, scal)
    }
    muMat <- modelMeans(Models, doses, FALSE, off, scal)
  }
  nD <- length(doses)
  Ntype <- match.arg(typeN)
  if(!is.null(alRatio)){
    if(length(alRatio)!= nrow(muMat)){
      stop("alRatio needs to be of same length as number of dose levels")
    }
    if(any(alRatio <= 0)){
       stop("all entries of alRatio need to be positive")
    }
    if (Ntype == "arm") {
       alRatio <- alRatio/min(alRatio)
    } else {
       alRatio <- alRatio/sum(alRatio)
    }
  } else {
    alRatio <- 1
  }
  contMat <- modContr(muMat, upperN * alRatio)       # no need to change
  summPow <-
    function(n, nD, contMat, muMat, sigma, alpha, sumFct, control, twoSide, ...)
    {
      ## returns "combined" power for a specified sample size n
      pow <- numeric(ncol(muMat))
      cVal <- critVal(contMat, n, alpha, control, twoSide)
      for(i in 1:ncol(muMat)) {
        pow[i] <- powCalc(contMat, n, cVal = cVal, mu = muMat[,i],
                          sigma = sigma, twoSide = twoSide, 
                          control = control)   
      }
      res <- sumFct(pow, ...)
      names(pow) <- dimnames(contMat)[[2]]
      attr(res, "power under models") <- round(pow, 4)
      res
    }
  ## need to ensure that limits on n bracket power
  upperPow <-
    summPow(round(upperN * alRatio), nD, contMat, muMat,
            sigma, alpha, sumFct, control, twoSide, ...)  
  if(upperPow < power){
    if(Ntype == "arm") sent <- "sample size in smallest group"
    else sent <- "total sample size"
    warning("upper limit for ", sent," is raised")
  }
  while (upperPow < power) {
    upperN <- 2 * upperN
    upperPow <-
      summPow(round(upperN * alRatio), nD, contMat, muMat,
              sigma, alpha, sumFct, control, twoSide, ...) 
  }
  lowerPow <-
    summPow(round(lowerN * alRatio), nD, contMat, muMat, sigma,
            alpha, sumFct, control, twoSide, ...)
  if(lowerPow > power){
    if(Ntype == "arm") sent <- "sample size in smallest group"
    else sent <- "total sample size"
    warning("lower limit for ", sent," is decreased")
  }
  while (lowerPow > power) {
    lowerN <- round(lowerN/2)
    if (lowerN == 0) stop("cannot find lower limit on n")
    lowerPow <-
      summPow(round(lowerN * alRatio), nD, contMat, muMat, sigma,
              alpha, sumFct, control, twoSide, ...)
  }
  ## now use simple binary search
  if (verbose) {
    cat("Upper N:", upperN, "Upper power", round(upperPow, 3), "\n")
    cat("Lower N:", lowerN, "Lower power", round(lowerPow, 3), "\n\n")
  }
  currPow <- 0
  niter <- 0
  while (abs(currPow - power) > tol & (upperN > lowerN + 1)) {
    currN <- round((upperN + lowerN)/2)
    currPow <- 
      summPow(round(currN * alRatio), nD, contMat, muMat, sigma,
              alpha, sumFct, control, twoSide, ...)
    if (currPow > power) {
      upperN <- currN
    } else {
      lowerN <- currN
    }
    niter <- niter + 1
    if (verbose){
      cat("Iter: ", niter, ", N = ", currN, ", power = ", round(currPow, 3),
          "\n", sep = "")
    }
  }
  while (currPow < power) { # step up until currPow >= power
    currN <- currN + 1
    currPow <-
      summPow(round(currN * alRatio), nD, contMat, muMat, sigma,
              alpha, sumFct, control, twoSide, ...)
  }
  if(length(alRatio) > 1) {   # unequal allocation
     if (Ntype == "arm") {
        res <- list(samp.size = round(currN * alRatio),
                    approx.power = round(currPow, 4))
     } else {
        res <- list(samp.size = round(currN * alRatio),
                    approx.power = round(currPow, 4))
     }
  } else {
    res <- list(samp.size = currN, approx.power = round(currPow, 4))
  }
  attr(res, "alRatio") <- round(alRatio/min(alRatio), 4)
  attr(res, "power") <- power
  attr(res, "twoSide") <- twoSide
  attr(res, "alpha") <- alpha
  attr(res, "sumFct") <- deparse(substitute(sumFct))
  oldClass(res) <- "sampSize"
  res
}

print.sampSize <- function(x, ...){
  cat("MCPMod sampSize\n")
  cat("\n")
  cat("Input parameters:\n")
  cat(" Summary Function:", attr(x, "sumFct"), "\n")
  cat(" Desired combined power value:", attr(x, "power"),"\n")
  if(attr(x, "twoSide")) side <- "two-sided"
  else side <- "one-sided"
  cat(" Level of significance:", paste(attr(x, "alpha"), " (", side, ")", sep=""), "\n")
  alRatio <- attr(x, "alRatio")
  if(length(alRatio) == 1){
    alRatio <- "balanced"
    sampMsg <- "Sample size per group:"
  } else {
    sampMsg <- "Sample sizes:"
  }
  cat(" Allocations:", alRatio, "\n\n")
  cat(paste(sampMsg), paste(x$samp.size), "\n")
  cat("\n")
  cat("Associated", attr(x, "sumFct"),"power:", paste(x$approx.power),"\n")
  cat("Power under models:","\n") 
  print(attr(x$approx.power, "power under models"))
  cat("\n")
}    

### example call
### doses <- c(0,10,25,50,100,150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential=c(85),
###                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=T, nrow=2)) 
### sampSize(models, doses, base = 0, maxEff = 0.4, sigma = 1, 
###          alpha = 0.05, upperN = 100, scal = 200)
### sampSize(models, doses, base = 0, maxEff = 0.4, sigma = 1,
###        alpha = 0.05, upperN = 400, scal = 200, alRatio=c(2,1,1,1,1,1)) 
### dvec <- c(0, 10, 50, 100) 
### mu1 <- c(1, 2, 2, 2)
### mu2 <- c(1, 1, 2, 2)
### mu3 <- c(1, 1, 1, 2)
### mMat <- cbind(mu1, mu2, mu3)
### dimnames(mMat)[[1]] <- dvec 
### sampSize(muMat = mMat, doses = dvec, sigma = 1,
###        alpha = 0.05, upperN = 400, alRatio=c(2,1,1,1))


LP <- function(models, model, type = c("both", "LP1", "LP2"), paramRange, doses,
               base, maxEff, sigma, n,
               len = c(10, 1), nr = 1, alpha = 0.025, twoSide = FALSE, 
               off = 0.1*max(doses), scal = 1.2*max(doses), control = mvtnorm.control()){
  ## function to calculate loss in power associated with
  ## prior parameter mis-specification
  ## models, doses, base, maxEff, sigma as above
  ## model - model for which LP should be investigated
  ## type - one of LP1 or LP2 (see JBS paper for definition)
  ## paramRange - vector of lower and upper limit for the prior parameter 
  ##       for models with two prior parameters matrix with the ranges
  ##       for each prior parameter in the rows
  ## n - number (in case of balanced sample size) 
  ##      or vector of group sample sizes
  ## len - number of points between min(paramRange) and max(paramRange) 
  ##       on which LP is calculated has to be of length 2 in case of 
  ##       models with 2 prior parameter
  ## nr - in case there is more than one model from one model class 
  ##   in the candidate set (e.g. two emax models), nr gives the number of the model 
  ## alpha, twoSide, scal, off as above
  
  type <- match.arg(type)
  builtIn <- c("emax", "exponential", "quadratic",
               "betaMod", "logistic", "sigEmax")
  usedM <- names(models)
  allowedMod <- intersect(builtIn, usedM)
  if(!is.element(model, allowedMod)){
     msg <- paste("model must be one of", paste(allowedMod, collapse=", "), "\n")
     stop(msg)
  }
  if (inherits(models, "fullMod")) {
    doses <- attr(models, "doses")
    base <- attr(models, "base")
    maxEff <- attr(models, "maxEff")
    off <- attr(models, "off")
    scal <- attr(models, "scal")
    Models <- models
  } else {
    Models <- fullMod(models, doses, base, maxEff, off, scal)
  }
  nD <- length(doses)
  if(length(n)==1) n <- rep(n,nD)
  nDF <- sum(n) - nD
  # calculate contrasts and critical value
  pM <- planMM(Models, doses, n, off, scal, FALSE, alpha, twoSide, control)
  contMat <- pM$contMat
  critV <- pM$critVal
  corMat <- pM$corMat
  # calculate mean vector under 'model'
  # determine, whether there is more than one prior estimate for model "model"
  if(nrInd <- length(Models[[model]]) > 4) mod <- list(Models[[model]][nr,])
  else mod <- list(Models[[model]])
  names(mod) <- model   
  mu <- modelMeans(mod, doses, std=FALSE, off, scal)
  # adjust arguments
  twoPars <- is.element(model, c("logistic", "betaMod", "sigEmax"))
  if(!twoPars){ # put paramRange in the same format
    paramRange <- matrix(c(paramRange, 0, 0), nrow=2, byrow=TRUE)
    len <- c(len,1)
  }
  if(type[1] != "LP1"){ # "LP2" or "both"
    ModelsIt <- Models
  }
  if(type[1] != "LP2"){ # "LP1" or "both"
    # calculate nominal power
    nomPower <- powCalc(contMat, n, alpha, mu = mu, sigma = sigma, 
                  cVal = critV, corMat = corMat, 
                  twoSide = twoSide, control = control)[1]
  }
  actPower <- potPower <- numeric(prod(len))
  j <- 1
  # parameter ranges
  sq1 <- seq(min(paramRange[1,]),max(paramRange[1,]),length=len[1])
  sq2 <- seq(min(paramRange[2,]),max(paramRange[2,]),length=len[2])
  
  for(h in sq2) {
    for(i in sq1) {
      if(twoPars) mod <- list(c(i,h))
      else mod <- list(i)
      names(mod) <- model
      mod <- fullMod(mod, doses, base, maxEff, off, scal)
      # mean vector under model with prior parameter i,h
      mu <- modelMeans(mod, doses, std=FALSE, off, scal)
      actPower[j] <- powCalc(contMat, n, alpha, mu = mu, sigma = sigma,
                       cVal = critV, corMat = corMat,
                       twoSide = twoSide, control = control)
      if(type[1] == "LP2" | type[1] == "both"){  # recalculate opt. cont. and critical value
        if(twoPars){          # using i,h as prior parameters
          if(nrInd) ModelsIt[[model]][nr,3:4] <- c(i,h)
          else ModelsIt[[model]][3:4] <- c(i,h)
        } else {
          if(nrInd) ModelsIt[[model]][nr,3] <- i
          else ModelsIt[[model]][3] <- i         
        }
        pMIt <- planMM(ModelsIt, doses, n, off, scal, FALSE, alpha, twoSide)
        potPower[j] <- powCalc(pMIt$contMat, n, alpha, mu = mu, sigma = sigma,
                        cVal = pMIt$critVal, corMat = pMIt$corMat, 
                        twoSide = twoSide, control = control)
      }
      j <- j + 1
    }
  }
  # determine parametername for output
  par <- switch(model,
          logistic = c("ED50", "delta"),
          betaMod = c("delta1", "delta2"),
          sigEmax = c("ED50", "h"),
          emax = "ED50",
          "delta")
  p1 <- rep(sq1,len[2])
  if(twoPars){
    p2 <- rep(rep(round(sq2,2),length=len[2]), each=len[1]) # values of 2nd parameter
    if(nrInd) used <- Models[[model]][nr,3:4] # determine used value
    else  used <- Models[[model]][3:4]
  } else {
    if(nrInd) used <- Models[[model]][nr,3]
    else  used <- Models[[model]][3]
  }
  res <- cbind(actPower) # matrix for output
  if(type[1] == "LP1"){
    LP1 <- nomPower - actPower
    res <- cbind(res, LP1)
  } else if(type[1] == "LP2"){
    res <- cbind(potPower, res)
    LP2 <- potPower - actPower
    res <- cbind(res, LP2)    
  } else if(type[1] == "both"){
    res <- cbind(potPower, res)
    LP1 <- nomPower - actPower
    LP2 <- potPower - actPower
    res <- cbind(res, LP1, LP2)      
  }
  if(twoPars){  # names for prior parameters
    res <- cbind(p1, p2, res)
    dimnames(res)[[2]][1:2] <- par
  } else {
    res <- cbind(p1, res)
    dimnames(res)[[2]][1] <- par    
  }
  attr(res, "model") <- model
  attr(res, "len") <- len
  attr(res, "used") <- used
  attr(res, "sampSizes") <- n
  attr(res, "type") <- type[1]
  attr(res, "twoPars") <- twoPars
  if(type[1] != "LP2") attr(res, "nomPower") <- nomPower
  oldClass(res) <- "LP"
  res
}
  
print.LP <- function(x, digits = 3, ...){
  cat("MCPMod LP","\n")
  cat("\n")
  cat("Model:", paste(attr(x, "model")),"\n")
  cat("Used Prior Parameter:", paste(attr(x, "used")), "\n")
  cat("Sample Sizes:", paste(attr(x, "sampSizes")),"\n")
  type <- attr(x, "type")
  if(type=="LP1" | type=="both"){
    cat("Nominal Power:", paste(round(attr(x, "nomPower"), digits)),"\n")
  }
  cat("\n")
  cat("Calculated Power Values:", "\n")
  print(round(data.frame(unclass(x)), digits = digits))
}

plot.LP <- function(x, line = TRUE, type = NULL, spldf = 5, ...){
  ## line - logical determining whether smooth line
  ##        should be fit to power values
  ## type - one of "LP1","LP2","both" or NULL
  ##        if NULL the type of the LP call is used

  if(line){      
    val <- 1
    len <- attr(x, "len")
    if(len[1] < 4)
      stop("at least 4 points needed to obtain smooth curve. Try: line = FALSE")      
    if (!(spldf > 1 && spldf <= dim(x)[1]))
      warning("you must supply 1 < spldf <= len")
  }
  twoPars <- attr(x, "twoPars")
  model <- attr(x, "model")
  used <- attr(x, "used")
  par <- dimnames(x)[[2]][1:(1+twoPars)]
  par1 <- x[,1]
  if(is.na(par[2])){
    par2 <- rep(1, dim(x)[1])
  } else {
    par2 <- x[,2]
  }
  if(is.null(type)){  # check if type is given in LP matrix
    type <- attr(x, "type") 
  } else {
    if(attr(x, "type")!="both" & attr(x, "type")!=type){
      stop("invalid type selected")
    }
  }
  # matrix for plotting data
  if(type == "both"){
    type <- c("LP1", "LP2")
    grp <- rep(type, each = dim(x)[1])
    pData <- data.frame(LP = c(x[, type]), par1 = rep(par1, 2), 
                 par2 = rep(paste(paste(par[2],"="), par2), 2), 
                 group = grp)
  } else {
    grp <- rep(type, dim(x)[1])
    pData <- data.frame(LP = x[, type], par1 = par1, 
                 par2 = paste(paste(par[2],"="), par2), 
                 group = grp)
  }
  main <- list(paste("Model:", model, ", Used value:", paste(used, collapse=" ")))
  L <- length(unique(grp))
  spL <- trellis.par.get("superpose.symbol")
  col <- spL$col[1:L] # use the dot color for the lines
  if(L == 2){ # key only needed when both LP1 and LP2 are to be plotted
    keyList <- list(points = list(col = col, pch = 1), transparent = TRUE,      
        text = list(paste(type), cex = 0.9), columns =  2)
    ylab <- NULL
  } else {
    keyList <- NULL
    ylab <- type
  }
  xyplot(LP~par1|par2, data=pData, 
         panel.data = list(uv = used[1], line = line, col = col, df = spldf),
         groups = grp, 
         panel = function(x, y, subscripts, groups,..., panel.data) {
           panel.abline(h = 0, lty = 2)
           panel.abline(v = panel.data[["uv"]], lty = 2)
           s <- seq(min(x), max(x), length = 101)
           panel.superpose(x, y, type = "p", lwd = 1.3, groups = groups, subscripts = subscripts)

           j <- 1
           if(panel.data[["line"]]){   
             for(i in unique(groups[subscripts])){   # add smooth line
               ind <- which(groups[subscripts]==i)
               p <- predict(smooth.spline(x[ind], y[ind], df = panel.data[["df"]]), s)
               panel.xyplot(p$x, p$y , type = "l", col = panel.data[["col"]][j])
               j <- j + 1
             } 
           }
         }, 
        main = main, xlab=par[1], strip=twoPars, ylab = ylab,
        key = keyList)
}


### example call
### doses <- c(0,10,25,50,100,150)
### models <- list(linear=NULL, emax=c(25), 
###            logistic=c(50,10.88111), exponential=c(85), 
###            betaMod=matrix(c(0.33,2.31,1.39,1.39),byrow=T,nrow=2))
### 
### # Examples from JBS paper, p.654
### LPobj <- LP(models, model = "emax", type = "both", paramRange = c(10,70),
###            doses = doses, base = 0, maxEff = 0.4, sigma = 1, n = 60,
###            alpha = 0.05, len = 15, scal = 200)
### print(LPobj)
### plot(LPobj)
###
### LPobj <- LP(models, model = "emax", type = "LP1", paramRange = c(10,70),
###           doses = doses, base = 0, maxEff = 0.4, sigma = 1, n = 60,
###           alpha = 0.05, len = 15, scal = 200)
### print(LPobj)
### plot(LPobj)
###
### # for exponential model with fullMod:
### fMod <- fullMod(models, doses, 0, 0.4, scal=200)
### LP(fMod, "exponential", "LP2", c(50, 120), sigma = 1, alpha = 0.05,
###    len = 20, n = 60)
### LPobj <- LP(fMod, "exponential", "both", c(50, 120), sigma = 1, alpha = 0.05,
###             len = 20, n = 60)
### plot(LPobj)
###
### Examples for models with two prior parameters
### LPobj <- LP(models, model = "betaMod", type = "both",
###            paramRange = matrix(c(0.3,1.9,0.4,2.5),nrow=2),
###            doses, base = 0, maxEff = 0.4, sigma = 1, n = 60, alpha = 0.05,
###            len = c(10,4), scal = 200)
### print(LPobj)
### plot(LPobj)
###
### LP(models, model = "logistic", type = "LP2",
###     paramRange = matrix(c(40,8,60,13),nrow=2),
###     doses, base = 0, maxEff = 0.4, sigma = 1, n = 60, 
###     alpha = 0.05, len = c(10,4), scal=200)
