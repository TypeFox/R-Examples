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

## functions to analyze dose finding data according to MCP-Mod 

## t-stats for MCP step
getTstat <-
  function(data, n, contMat, resp = "resp", dose = "dose")
{
  if (any(is.na(match(c(resp, dose), names(data))))) {
    stop(resp," and/or ", dose, " not found in data")
  }
  dose <- data[, dose]
  resp <- data[, resp]
  ## means per dose
  mn <- tapply(resp, dose, mean)

  ## pooled standard deviation
  sdv <- tapply(resp, dose, sd)
  if (length(n) == 1) n <- rep(n, length(sdv))
  # remove NAs for dose groups with only 1 obs.
  if(any(n==1)){
    sdv[n==1] <- 0
  }
  n1 <- n - 1
  sdv <- sqrt(sum(n1 * sdv^2)/sum(n1))

  ## contrasts
  ct <- as.vector(mn %*% contMat)
  den <- sdv * sqrt(colSums((contMat^2)/n))

  ## max t-stat
  val <- ct/den
  names(val) <- dimnames(contMat)[[2]]
  val
}

# function to calculate p-values
pValues <-
  function(cMat, n, alpha = 0.025, Tstats, control = mvtnorm.control(),
           twoSide = FALSE, corMat = NULL, ...)
{
  ## function to calculate p-values
  ##
  ## cMat - model contrast matrix nDose x nMod
  ## n - scalar or vector with sample size(s)
  ## alpha - significance level
  ## control - see mvtnorm.control function
  ## twoSide - logical variable indicating if one or two-sided test
  ## corMat - correlation matrix

  nD <- nrow(cMat)
  nMod <- ncol(cMat)
  if(length(Tstats) != nMod){
    stop("Tstats needs to have length equal to the number of models")
  }
  if (length(n) == 1) {
    nDF <- nD * (n - 1)
  } else {
    if (length(n) != nD) {
      stop("'n' must have length as number of doses")
    }
    nDF <- sum(n) - nD
  }
  if(is.null(corMat)){
      corMat <- t(cMat)%*%(cMat/n)
      den  <- sqrt(crossprod(t(colSums(cMat^2/n))))
      corMat <- corMat / den
  }
  ctrl <- mvtnorm.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  if(twoSide){
    lower <- matrix(rep(-Tstats, each = nMod), nrow = nMod)
  }
  else lower <- matrix(rep(-Inf, nMod^2), nrow = nMod)
  upper <- matrix(rep(Tstats, each = nMod), nrow = nMod)
  pVals <- numeric(nMod)
  for(i in 1:nMod){
    pVals[i] <-   1 - pmvt(lower[,i], upper[,i], df = nDF, corr = corMat, 
                   algorithm = ctrl, ...)
  }
  pVals
}

### functions to fit models, and functions to get starting values
### for the nls function

### Initial estimates for built-in models
### two auxiliary functions to get starting values for the nls algorithm
getInitLRasymp <-
  ## left and right asymptotes
  function(meanVal, dlt = 0.1)
{
  rg <- range(meanVal)
  dlt <- dlt * diff(rg)
  c(e0 = rg[1] - dlt, eMax = rg[2] + dlt)
}

getInitP <-
  ## dose that gives certain percentage of max effect
  function(meanVal, eVec = getInitLRasymp(meanVal), p = 0.5, dose, ve=FALSE)
{
  e0   <- eVec[1]
  eMax <- eVec[2]
  targ <- e0 * (1 - p) + eMax * p
  ind  <- meanVal <= targ
  ind1 <- !ind
  if (!any(ind) || all(ind)) stop("cannot find initial values")
  ind  <- ind & (dose <= max(dose[ind1]))
  d1   <- max(dose[ind]); p1 <- meanVal[dose == d1]
  ind1 <- ind1 & dose > d1
  d2   <- min(dose[ind1]); p2 <- meanVal[dose == d2]
  res <- d1 + (d2 - d1) * (targ - p1)/(p2 - p1)
  names(res) <- NULL
  res
}

getInitBeta <-
  function(m, scal, dose)
{
  ## calculates starting values for the beta model
  ## for the nls algorithm
  ## dose and m assumed to be ordered 
  ## in the same order (according to dose)
   dmax <- dose[which(m==max(m))]
   e0 <- m[dose==0]
   emax <- max(m) - e0
   ds <- dose[order(m)[length(m)-1]] # select as 2nd dose the one
                                     # with 2nd highest resp.
   z <- dmax/(scal-dmax)
   f <- (m[dose==ds] - e0)/emax
   beta <- try(log(f)/(log(ds^z*(scal-ds)) - log(dmax^z*(scal-dmax))))
   alpha <- z*beta
   if(is.na(alpha)) alpha <- beta <- 1 # may occur if f < 0; (1,1) as
                                       # initial values seems always to
                                       # work quite well 
   res <- c(alpha, beta)
   names(res) <- NULL
   res
}

getInit <-
  ## Initial estimates for built-in models, if needed
  ## non-linear fitting is done with p-linear, but Gauss-Newton used
  ## later due to some problems with se.fit of the predictions of the
  ## p-linear fit 
  function(data, model = c("emax", "exponential", "logistic", "betaMod",
                   "sigEmax"), scal)
{
  model <- match.arg(model)             # only built-in allowed
  meanVal <- data$respM
  dose <- data$doseM
  
  switch(model,
         emax = {
           c(ed50 = getInitP(meanVal, dose = dose))
         },
         exponential = {
           e0 <- getInitLRasymp(meanVal)[1]
           meanVal <- meanVal - e0
           aux <- coef(lm(log(meanVal) ~ dose, na.action = na.omit))
           names(aux) <- NULL
           c(delta = 1/aux[2])
         },
         logistic = {
           ed50 <- getInitP(meanVal, dose = dose)
           delta <- getInitP(meanVal, p = 0.75, dose = dose) - ed50
           c(ed50 = ed50, delta = delta)
         },
         betaMod = {
           init <- getInitBeta(meanVal, scal, dose)
           c(delta1 = init[1], delta2 = init[2])
         },
         sigEmax = {
           ed50 <- getInitP(meanVal, dose = dose)
           ed10 <- getInitP(meanVal, p = 0.1, dose = dose)
           ed90 <- getInitP(meanVal, p = 0.9, dose = dose)
           h <- 1.91/log10(ed90/ed10) # From Ch. 14, Ting (2006),
                                      # Dose finding in drug development
           c(ed50 = ed50, h = h)
         })
}

fitModel <-
  ## fits dose-response model from list of built-in models,
  ## or according to model provided by user
  function(data, model, resp = "resp", dose = "dose", start,
           control = NULL, off = 1, scal = 1, uModPars = NULL,
           addArgs = NULL, na.action = na.fail)
{
  ## data - data.frame with information on dose and response
  ##          and possibly other covariates used in model
  ## model - string with name of model function to be used
  ## resp, dose - characters with names of columns in data
  ##              corresponding to response and dose
  ## control - optional list with control parameters for fit
  ## off - offset for dealing with dose = 0 in lin-log model
  ## uModPars - optional character vector with names/expressions
  ##            of user-defined model parameters (names(start) used by
  ##            default) 
  ## addArgs - optional character vector with names of additional
  ##           arguments (variables) to user-defined model
  ## na.action - function specifying how to handle NAs in data

  ## ensuring resp and dose are defined in data
  data$dose <- data[, dose]         
  data$resp <- data[, resp]

  ## reacting to possible NAs
  ## it would be better to only take into account NAs in
  ## varibles actually used in the model.
  data <- na.action(data)

  ## checking if built-in model and assigning it number, if so
  builtIn <- c("linlog", "linear", "quadratic", "emax",
        "exponential","logistic", "betaMod", "sigEmax")
  modelNum <- match(model, builtIn)

  ## under the built-in models with no covariates, the mean
  ## modeling can be used in all cases
  respM <- tapply(data$resp, data$dose, mean) # fit models based on means
  doseM <- sort(unique(data$dose))
  n <- as.vector(table(data$dose))
  dataM <- data.frame(doseM = doseM, respM = respM, n = n)
  vars <- tapply(data$resp, data$dose, var)
  # allow for dose groups with only one patient
  vars[n==1] <- 0 # replace NAs with 0
  S2 <- sum((n - 1) * vars)
  
  if (!is.na(modelNum)) {               # built-in model
    if (is.element(modelNum, 1:3)) {    # linear models
      if (modelNum == 1) {
        dataM$off <- rep(off, nrow(dataM))  
        ## linear-log model
        fm <- lm(respM ~ I(log(doseM + off)), dataM, qr=TRUE, weights = n)
        base <- coef(fm)[1]+coef(fm)[2]*log(off)
      }
      if (modelNum == 2) {
        ## linear model
        fm <- lm(respM ~ doseM, dataM, qr=TRUE, weights = n)
        base <- coef(fm)[1]
      }
      if (modelNum == 3) {
        ## quadratic model
        fm <- lm(respM ~ doseM + I(doseM^2), dataM, qr=TRUE, weights = n)
        base <- coef(fm)[1]               # baseline
      }
    } else {                            # nonlinear built-in model
      if (is.null(start)) {             # need to derive initial values
        start <- getInit(dataM, model, scal)
      }
      if (modelNum == 4) {      # Emax model
        names(start) <- NULL
        start <- c(led50 = log(start))
        fm <- try(nls(respM ~ cbind(1, emax(doseM, 0, 1, exp(led50))), 
                     start = start, data = dataM, model = TRUE,
                     algorithm = "plinear", control = control,
                     weights = n), silent = TRUE) 
        if (!inherits(fm, "nls")) {
          fm <- NA
        } else {
          base <- coef(fm)[2]
        }
      }
      if (modelNum == 5) {               # exponential model
        names(start) <- NULL
        start <- c(ldelta = log(start))
        fm <-
           try(nls(respM ~ cbind(1, exponential(doseM, 0, 1, exp(ldelta))),
                   start = start, weights = n, data = dataM, model = TRUE,
                   control = control, algorithm = "plinear"), silent = TRUE)
        if (!inherits(fm, "nls")) {
          fm <- NA
        } else {
          base <- coef(fm)[2]
        }
      }
      if (modelNum == 6) {               # logistic model
        start <- c(log(start["ed50"]), start["delta"])
        names(start) <- c("led50", "delta")
        fm <- try(nls(respM ~ cbind(1, logistic(doseM, 0, 1, exp(led50),
                                                delta)),
                      start = start, algorithm = "plinear", data = dataM,
                      model = TRUE, control = control, weights = n),
                      silent = TRUE)
        if (!inherits(fm, "nls")) {
          fm <- NA
        } else {
          base <- predict(fm, newdata = list(doseM = 0))
        }
      }
      if (modelNum == 7) {              # beta model
        dataM$scal <- rep(scal, nrow(dataM))  
        fm <-
          try(nls(respM ~ cbind(1, betaMod(doseM, 0, 1, delta1, delta2, scal)),
                  start = start, weights = n, data = dataM, model = TRUE,
                  algorithm = "plinear", control = control), silent = TRUE)
        if (!inherits(fm, "nls")) {
          fm <- NA
        } else {
          base <- coef(fm)[3]
        }
      }
      if (modelNum == 8) {              # sigEmax model
        start <- c(log(start["ed50"]), start["h"])
        names(start) <- c("led50", "h")
        fm <-
          try(nls(respM ~ cbind(1, sigEmax(doseM, 0, 1, exp(led50), h)),
                  start = start, model = TRUE, 
                  algorithm = "plinear", control = control, weights = n),
                  silent = TRUE)  
        if (!inherits(fm, "nls")) {
          fm <- NA
        } else {
          base <- coef(fm)[3]
        }
      }
    }
    if (length(fm) > 1) ResidSS <- S2 + deviance(fm)
  } else {                              # user defined model
    ## first check for starting estimates, parameter names, etc
    if (is.null(start)) {
      stop("must provide stating estimates for user-defined model")
    }
    namStart <- names(start)
    if (is.null(namStart)) {
      stop("'start' must have names for user-defined models")
    }
    if (is.null(uModPars)) uModPars <- namStart  # default parameters
    ## create the model call
    modForm <- paste("resp ~ ", model, "(dose,",
                     paste(uModPars, collapse = ","), sep = "")
    if (!is.null(addArgs)) {
      modForm <- paste(modForm, ",", paste(addArgs, collapse = ","))
    }
    modForm <- paste(modForm, ")")
    modForm <- eval(parse(text = modForm))
    ## will assume non-linear model
    fm <- try(do.call("nls", list(modForm, data, start, control)))
    if (!inherits(fm, "nls")) {
      fm <- NA
    } else {
      ResidSS <- deviance(fm)
      dataM <- NULL
      ## following will not work when covariates are used in nls model
      base <- predict(fm, newdata = list(dose = 0))
    }
  }
  if (length(fm) > 1) {
    res <- list(fm = fm, base = base)
    ## saving information for other methods
    attr(res, "ResidSS") <- ResidSS
    attr(res, "model") <- model
    attr(res, "scal") <- scal
    attr(res, "off") <- off
    attr(res, "dataM") <- dataM
  } else {
    res <- NA
  }
  oldClass(res) <- "fitMCPMod"          # allow use of methods later
  res
}

# function to return analyt. gradient for built-in or
# user defined models
getGrad <-
  function(obj, dose, uGrad = NULL)
  ## obj - an object inheriting from class fitMCPMod
  ## dose - vector with doses at which to calculate the  gradient
  ## uGrad - optionally, a character string with the name of a
  ##         user-defined gradient function
{
  if(!inherits(obj, "fitMCPMod")) {
    stop("obj must inherit from class fitMCPMod")
  }
  model <- attr(obj, "model")
  fm <- obj$fm
  if (!inherits(fm, "nls")) {
    ## linear models are left out
    stop("fitted object must inherit from class nls")
  }
  ## only used internally in getGrad
  lg2 <- function(x) ifelse(x == 0, 0, log(x))

  cf <- coef(obj)                       # estimated parameters
  res <-
    switch(model,
           "emax" = {
             eMax <- cf[2]; ed50 <- cf[3]
             cbind(1, dose/(ed50 + dose), -eMax*dose/(dose + ed50)^2)
           },
           "logistic" = {
             eMax <- cf[2]; ed50 <- cf[3]; delta <- cf[4]
             den <- 1 + exp((ed50 - dose)/delta)
             g1 <- -eMax*(den-1)/(delta*den^2)
             g2 <- eMax*(den-1)*(ed50-dose)/(delta^2*den^2)
             cbind(1, 1/den, g1, g2)
           },
           "sigEmax" = {
             eMax <- cf[2]; ed50 <- cf[3]; h <- cf[4]
             den <- (ed50^h + dose^h)
             g1 <- dose^h/den
             g2 <- -ed50^(h-1)*dose^h*h*eMax/den^2
             g3 <- eMax*dose^h*ed50^h*lg2(dose/ed50)/den^2
             cbind(1, g1, g2, g3)
           },
           "betaMod" = {
             scal <- attr(obj, "scal")
             dose <- dose/scal
             if (any(dose > 1)) {
               stop("doses cannot be larger than scal in betaModel")
             }
             delta1 <- cf[3]; delta2 <- cf[4]; eMax <- cf[2] 
             maxDens <- (delta1^delta1)*(delta2^delta2)/
               ((delta1 + delta2)^(delta1+delta2))
             g1 <- ((dose^delta1) * (1 - dose)^delta2)/maxDens
             g2 <- g1*eMax*(lg2(dose)+lg2(delta1+delta2)-lg2(delta1))
             g3 <- g1*eMax*(lg2(1-dose)+lg2(delta1+delta2)-lg2(delta2))
             cbind(1, g1, g2, g3)
           },
           "exponential" = {
             delta <- cf[3]
             e1 <- cf[2]
             cbind(1, exp(dose/delta) - 1, -exp(dose/delta)*dose*e1/delta^2 )
           },
           {
             ## user defined gradient function
             if(is.null(uGrad)) {
               stop("user-defined gradient needs to be specified")
             }
             do.call(uGrad, c(list(dose), cf))
           })
  res
}

# function to return coefficients from fit
coef.fitMCPMod <-
  function(object, ...)
{
  ## object - object inheriting from class MCPMod
  fm <- object$fm
  cf <- coef(fm)
  if(inherits(fm, "lm")) return(cf)     
  temp <- names(cf) # save names for user-models
  names(cf) <- NULL
  model <- attr(object, "model")
  cf <- switch(model,
          "emax" = c(e0 = cf[2], eMax = cf[3], ed50 = exp(cf[1])),
          "logistic" = c(e0=cf[3], eMax=cf[4], ed50 = exp(cf[1]),
            delta = cf[2]),
          "exponential" = c(e0 = cf[2], e1 = cf[3], delta = exp(cf[1])),
          "betaMod" = c(e0=cf[3], eMax=cf[4], delta1=cf[1], delta2=cf[2]),
          "sigEmax" = c(e0=cf[3], eMax=cf[4], ed50=exp(cf[1]), h=cf[2]),
          {
          names(cf) <- temp
          cf
          }
        )
  cf
}

# predict method
predict.fitMCPMod <-
  function(object, doseSeq, se.fit = TRUE, uGrad = NULL, ...)
{
  # object - either a lm object or a nls object in the latter case
  #     the code distinguishes between user-defined models
  #     and built-in models
  # doseSeq - sequence for predictions
  # se.fit - logical indicating whether sd of the mean should be
  #          calculated
  model <- attr(object, "model")
  scal <- attr(object, "scal")
  off <- attr(object, "off")
  fm <- object$fm
  if(inherits(fm, "lm")){ # in this case use the implemented predict method
    newdata <- data.frame(doseM = doseSeq, off = rep(off, length(doseSeq)))
    res <- predict(fm, newdata, se.fit = se.fit)
    ## Need to correct se.fit, if requested
    if(se.fit) {
      df <- sum(fm$weights) - length(coef(fm))
      sigCorr <- sqrt(attr(object, "ResidSS")/df)
      sigOrig <- summary(fm)$sigma
      res$se.fit <- (sigCorr/sigOrig) * res$se.fit
      res$df <- df
      res$residual.scale <- sigCorr
    }
    return(res)
  }
  builtIn <- is.element(model, c("emax", "exponential", "logistic",
                              "betaMod", "sigEmax"))
  if(inherits(fm, "nls") & !builtIn){ # user defined model
    mn <- predict(fm, data.frame(dose = doseSeq))
    if(se.fit) {
      ## first, check if model includes a gradient attribute
      grd <- attr(mn, "gradient")
      if(is.null(grd)) {
        if(is.null(uGrad)) {
          stop("neither gradient attribute, nor uGrad defined")
        }
        grd <- getGrad(object, doseSeq, uGrad)
      }
      Rinv <- solve(fm$m$Rmat())
      sig <- summary(fm)$sigma
      seFit <- sig * sqrt(rowSums((grd %*% Rinv)^2))
      df <- df.residual(fm)
      res <- list(fit = mn, se.fit =  as.vector(seFit),
                  residual.scale = sig, df = df)
      return(res)
    } else {
      return(mn)
    }
  } else { # built in model
    dd <- data.frame(doseM = doseSeq) # dose levels to be predicted
    cf <- coef(object)  # extract coefficients
    if(model == "betaMod"){ 
      cf <- c(cf, scal) # add scale parameter for betaModel 
      dd$scal <- scal
    }
    mn <- predict(fm, dd) # predict mean response
    if(!se.fit){
      return(mn)
    } 
    RSS <- attr(object, "ResidSS")   # get residual sum of squares
    doseV <- attr(object, "dataM")$doseM # recover dose-vector
    n <- fm$weights
    df <- sum(n) - length(cf) 
    sig <- sqrt(RSS/df) # residual variance estimate
    ## Gradient matrix: see Bates/Watts p. 58/59
    V <- getGrad(object, doseV)
    if(all(!is.infinite(V)) & all(!is.nan(V))){
      ## checking for infinities (can happen because the analytic
      ## gradient is not used for fitting)
      V <- t(crossprod(V, sqrt(diag(n))))  # weights for unequal sample sizes
      R <- qr.R(qr(V))
      Rinv <- try(solve(R))
      if (!inherits(Rinv, "matrix")){
        ## sometimes R is singular (typically for very unequally
        ## spaced dose designs)
        warning("Cannot cannot calculate standard deviation for ",
                    model, " model.\n") 
        seFit <- rep(NA, length(doseSeq))            
      } else {
        v <- getGrad(object, doseSeq) # calc. gradient
        seFit <- sig * sqrt(rowSums((v%*%Rinv)^2))     
      }
    } else {
      warning("Cannot cannot calculate standard deviation for ",
                    model, " model.\n") 
      seFit <- rep(NA, length(doseSeq))       
    }
    res <- list(fit = mn, se.fit =  as.vector(seFit),
                residual.scale = sig, df = df)
  }
  res
}

# Calculate information criteria for nls and lm objects
AIC.fitMCPMod <-
  function(object, ..., k = 2)
{
  # object - an object of class fitMCPMod
  # k - penalty term
  fm <- object$fm
  if (is.null(fm$weights)) {
    ## user defined model or non-averaged modeling
    return(AIC(fm, k = k))
  }
  RSS <- attr(object, "ResidSS")
  n <- sum(fm$weights)
  sig2 <- RSS/n
  logL <- -n/2*(log(2*pi) + 1 + log(sig2))
  -2*logL + k*(length(coef(object)) + 1) # "+ 1" because of sigma parameter
}

## model selection
modelSelect <-
  function(data, namSigMod, selMethod, pW, resp, dose, start, nlsControl,
           off, scal, uModPars, addArgs)
{ 
  fm <- NA
  warn <- NULL
  nSigMod <- length(namSigMod)
  if (selMethod == "maxT") {  # first maxT option
    i <- 1
    while((length(fm) == 1) && (i <= nSigMod)) {
      nam <- namSigMod[i]
      fm <- fitModel(data, nam, resp, dose, start[[nam]], nlsControl, 
                     off, scal, uModPars[[nam]], addArgs[[nam]])
      if(length(fm) == 1) # model didn't converge
        warning(nam, " model did not converge\n")
      i <- i + 1
    }
    if (length(fm) > 1) {  # model converged
      fm <- list(fit = list(fm), base = fm$base)
      attr(fm$fit, "model2") <- names(fm$fit) <- nam
    } else {
      fm <- list(fit = NA, base = NA) # no model converged 
    }
  } else {                    # AIC, BIC, aveAIC or aveBIC
    fm <- vector("list", nSigMod)
    crit <- fmBase <- rep(NA, nSigMod)
    if (selMethod == "AIC"| selMethod == "aveAIC") { 
      pen <- 2
    } else {
      pen <- log(dim(data)[[1]])
    }
    names(fm) <- names(crit) <- names(fmBase) <- namSigMod
    for(i in 1:nSigMod) {
      nam <- namSigMod[i]
      fitmod <- fitModel(data, nam, resp, dose, start[[nam]], nlsControl, 
                          off, scal, uModPars[[nam]], addArgs[[nam]])
      if(!is.list(fitmod)) { # model didn't converge
        fm[[i]] <-  NA
        warning(nam, " model did not converge\n")
      } else { # model converged
        fm[[i]] <- fitmod
        fmBase[i] <- fitmod$base
        crit[i] <- AIC(fitmod, k = pen)
      }
    }
    if (all(is.na(crit))) {
      fm <- NA
    } else { 
      if (selMethod == "AIC" | selMethod == "BIC") {
        model2 <- namSigMod[which(crit == min(crit, na.rm = TRUE))]
        fm <- list(fm[[model2]])
        fmBase <- fmBase[model2]
        attr(fm, "model2") <- names(fm) <- model2
        attr(fm, "IC") <- crit
      }
      else { # model averaging
        attr(fm, "model2") <- namSigMod[!is.na(fmBase)]
        attr(fm, "IC") <- crit
        crit <- crit[!is.na(crit)]
        ## subtract const from crit values to avoid numerically 0
        ## values (exp(-0.5*1500)=0!)
        const <- mean(crit)
        if(is.null(pW)){
          pW <- rep(1, length(crit)) # standard 'noninformative' prior
          names(pW) <- names(crit)
        } else {
        pW <- pW[names(crit)]
        if(any(is.na(pW))) stop("pW needs to be a named vector with names equal to the models in the candidate set")
        }
        attr(fm, "weights") <-
          pW*exp(-0.5*(crit-const))/sum(pW*exp(-0.5*(crit-const)))
        attr(fm, "pweights") <- pW
      }
    }
    fm <- list(fit = fm, base = fmBase)   
  }
  fm
}

## dose estimation
getDose <-
  function(dose, ind)
{
  aa <- !is.na(ind)
  if (!all(aa)) {
    ind <- ind[aa]; dose <- dose[aa]
  }
  if (any(ind)) min(dose[ind])
  else NA
}

getDoseEst <-
  function(fmb, clinRel, rangeDose, doseEst, selMethod, 
           dePar = 0.1, lenDose = 101, uGrad = NULL)
{
  ## equally spaced points within dose range for prediction
  doseSeq <- seq(rangeDose[1], rangeDose[2], len = lenDose)
  off <- attr(fmb, "off"); scal <- attr(fmb, "scal")
  newdata <- list(dose = doseSeq, off = rep(off, lenDose), scal = rep(scal, lenDose))
  ldePar  <- length(dePar)
  val <- double(ldePar)
  gma <- 1-dePar
  base <- fmb$base
  tDose <- matrix(ncol = length(base)-sum(is.na(base)), nrow = length(dePar))
  z <- 1
  for(m in which(!is.na(base))){ # only predict models that converged
    if(doseEst[1] != "ED"){ # MED estimate
      pred <- predict(fmb$fit[[m]], doseSeq, uGrad = uGrad)
      fo   <- data.frame(pred = pred$fit - base[m], se.fit = pred$se.fit)
      df   <- pred$df
      ## selects the desired dose according to MED method
      for(i in 1:ldePar) {
        crt  <- qt(gma[i], df)
        LL   <- fo$pred - crt * fo$se.fit
        UL   <- fo$pred + crt * fo$se.fit
        switch(doseEst,
             "MED1" = ind <- LL >= 0 & UL >= clinRel,
             "MED2" = ind <- LL >= 0 & fo$pred >= clinRel,
             "MED3" = ind <- LL >= clinRel
          )
        val[i] <- getDose(doseSeq, ind)
      }
    } else { # ED estimate
      pred <- predict(fmb$fit[[m]], doseSeq, se.fit = FALSE)
      pEff <- dePar*(max(pred) - base[m])
      for(i in 1:ldePar) {
        ind <- pred >= pEff[i] + base[m]
        val[i] <- getDose(doseSeq, ind)
      }
    }
    tDose[,z] <- val    
    z <- z + 1
  }
  if(is.element(selMethod, c("aveAIC", "aveBIC"))){
     doseAve <- as.vector(tDose%*%attr(fmb$fit, "weights"))
     if(doseEst[1] == "ED")
       names(doseAve) <- paste(doseEst, round(100 * dePar), "%", sep = "")
     else
       names(doseAve) <- paste(doseEst,"," , round(100 * (1-2*dePar)), "%",
                               sep = "")
     dimnames(tDose) <- list(names(doseAve), attr(fmb$fit, "model2"))
     attr(doseAve, "tdModels") <- tDose
     tDose <- doseAve
  } else {
     tDose <- as.vector(tDose)
     if(doseEst == "ED")                                                      
       names(tDose) <- paste(doseEst, round(100 * dePar), "%", sep = "")       
     else                                                                     
       names(tDose) <- paste(doseEst,"," , round(100 * (1-2*dePar)), "%",
                             sep = "")
  }
  tDose
} 

recovNames <-
  function(names)
{
## function to recover model names (in case of multiple models from one class)
## just for use in MCPMod function
## example: recovNames(c("emax1", "betaMod", "emax2", "logistic", "usermodel"))
## returns: c("emax","betaMod","logistic","usermodel")

   builtIn <- c("linlog", "linear", "quadratic", "emax",
        "exponential","logistic", "betaMod", "sigEmax")
   newnames <- character()
   i <- 1
   for(nam in names){
      pm <- pmatch(builtIn, nam)
      if(any(!is.na(pm))) {
        newnames[i] <- builtIn[which(!is.na(pm))]
        i <- i+1
      }
      if(all(is.na(pm))) {
        newnames[i] <- nam
        i <- i+1
      }
   }
   unique(newnames)
}

## function to get reasonable defaults for, off, scal and dePar
getDef <-
  function(off = NULL, scal = NULL, dePar = NULL, doses, doseEst)
{
  mD <- max(doses)
  if(is.null(dePar)){ # default for dePar depending on the dose estimate
    dePar <- ifelse(doseEst=="ED", 0.5, 0.1)
  }
  if(is.null(scal)){ # default for scal parameter
    scal <- 1.2*mD
  } else { # check if valid scal is provided
    if(scal < mD){
      stop("'scal' should be >= maximum dose")
    }
  }
  if(is.null(off)){ # default for off parameter
    off <- 0.1*mD
  }
  list(scal = scal, off = off, dePar = dePar)
}


### main function implementing methodology, combining several others
MCPMod <- 
  function(data, models = NULL, contMat = NULL, critV = NULL,
           resp = "resp", dose = "dose", off = NULL, scal = NULL,
           alpha = 0.025, twoSide = FALSE,
           selModel = c("maxT", "AIC", "BIC", "aveAIC", "aveBIC"),
           doseEst = c("MED2", "MED1", "MED3", "ED"),
           std = TRUE, start = NULL,
           uModPars = NULL, addArgs = NULL, 
           dePar = NULL, clinRel = NULL, lenDose = 101, pW = NULL, 
           control = list(maxiter = 100, tol = 1e-6, minFactor = 1/1024),
           signTtest = 1, pVal = FALSE, testOnly = FALSE,
           mvtcontrol = mvtnorm.control(), na.action = na.fail,
           uGrad = NULL)
{
  ## data - data.frame with, at least, columns "resp" and "dose"
  ## models - list with component names specifying models
  ##          and optionally elements being used to calculate
  ##          model contrasts; need to have named elements with
  ##          corresponding models and be consistent with contMat,
  ##          if that is non-null
  ## contMat - contrast matrix for candidate set and doses
  ## critV - critical value for MCP test
  ## alpha - significance level for calculating critVal (if = NULL)
  ## selModel - model selection criterion
  ## doseEst - type of dose estimator to be used
  ## dePar - gives "gamma" parameter for MED-type estimators or
  ##         the "p" in the EDp estimate
  ## pW - named vector of prior probs. for different models
  
  if (any(is.na(match(c(resp, dose), names(data))))) {
    stop(resp," and/or ", dose, " not found in data")   
  }
  data <- na.action(data)
  ind <- match(dose, names(data))
  data <- data[order(data[, ind]), ]
  ## target dose estimate type
  doseEst <- match.arg(doseEst)     
  n <- as.vector(table(data[, dose]))           # sample sizes per group
  doses <- sort(unique(data[, dose]))
  ## getting defaults which depend on the DRdata
  def <- getDef(off, scal, dePar, doses, doseEst) 
  scal <- def[[1]]; off <- def[[2]]; dePar <- def[[3]]
  if (is.null(contMat)) {
    ## need to calculate it
    mu <- modelMeans(models, doses, std, off, scal)
    contMat <- modContr(mu, n)
  }
  ## MCP test first
  tStat <- signTtest * getTstat(data, n, contMat, resp, dose)
  if(twoSide){
    tStat <- abs(tStat)
  }
  if (is.null(critV)) {
    if(pVal) {
      pVals <- pValues(contMat, n, alpha, tStat, mvtcontrol, twoSide)
    }
    critV <- critVal(contMat, n, alpha, mvtcontrol, twoSide = twoSide)  
    attr(critV, "Calc") <- TRUE # determines whether cVal was calculated
  } else { 
    pVal <- FALSE # pvals are not calculated if critV is supplied
    attr(critV, "Calc") <- FALSE
  } 
  indStat <- tStat > critV
  if (!any(indStat) | testOnly) {
    ## only mcp-test or no significant t-stats 
    result <- list(signf = any(indStat), model1 = NA, model2 = NA)
  } else {
    ## model selection method
    selMethod <- match.arg(selModel)
    control <- do.call("nls.control", control)
    ## at least one significant, select a model if possible
    namMod <- names(tStat)
    maxTstat <- max(tStat)
    model1 <- namMod[which(tStat == maxTstat)] # model with most sig contrast
    ## significant models, in descending order of tstat
    indSigMod <- 1:length(tStat)
    indSigMod <- indSigMod[rev(order(tStat))][1:sum(indStat)]
    namSigMod <- namMod[indSigMod]       # significant models
    namSigMod <- recovNames(namSigMod)   # remove model nrs.
                                         # (in case of multiple models)
    fmb <-
      modelSelect(data, namSigMod, selMethod, pW, resp, dose, start,
                  control, off, scal, uModPars, addArgs) 
    if (all(is.na(fmb$base))) { # none of sign. model converged
      result <- list(signf = TRUE, model1 = model1, model2 = NA)
    } else {
      ## dose estimation model(s) obtained
      ## move to final step, dose estimation
      ## dose range
      rgDose <- range(doses)
      tDose <- getDoseEst(fmb, clinRel, rgDose, doseEst, selMethod, 
                    dePar, lenDose, uGrad)
      result <- list(signf = TRUE, model1 = model1,
                     model2 = attr(fmb$fit, "model2"))
    }
  }
  ## add information to the object
  result$input <- list(models=models, resp=resp, dose=dose, off=off,
                       scal=scal, alpha=alpha, twoSide=twoSide,
                       selModel=selModel[1], doseEst=doseEst, 
                       std=std, dePar=dePar, uModPars=uModPars,
                       addArgs=addArgs, start = start, uGrad=uGrad,
                       clinRel=clinRel, lenDose=lenDose, signTtest=signTtest,
                       pVal=pVal, testOnly=testOnly)
  result$data <- data
  result$contMat <- contMat
  result$corMat <- t(contMat)%*%(contMat/n) /
    sqrt(crossprod(t(colSums(contMat^2/n))))
  result$cVal <- critV
  result$tStat <- tStat
  if(pVal) attr(result$tStat, "pVal") <- pVals
  if (!all(is.na(result$model2))){
    result$fm <- fmb$fit
    result$tdose <- tDose
  }
  oldClass(result) <- "MCPMod"
  result
}

# print method for MCPMod objects
print.MCPMod <-
  function(x, digits = 3,...)
{
  cat("MCPMod\n\n")
  if(x$input$twoSide) side <- "two-sided"
  else side <- "one-sided"
  if(!attr(x$cVal, "Calc")){
    cat(paste("PoC:",sep=""),
        paste(c("no", "yes")[x$signf+1]), "\n")    
  } else {
    cat(paste("PoC (alpha = ", x$input$alpha,", ", side, "):",sep=""),
      paste(c("no", "yes")[x$signf+1]), "\n")
  }
  if(x$signf & !x$input$testOnly){
    cat("Model with highest t-statistic:", paste(x$model1),"\n")
    if(!all(is.na(x$model2))){
      modSel <- is.element(x$input$selModel, c("AIC", "BIC", "maxT"))
      mo <- ifelse(modSel, "Model", "Models")
      cat(paste(mo, "used for dose estimation:"), paste(x$model2),"\n")
      cat("Dose estimate:","\n")
      attr(x$tdose, "tdModels") <- NULL 
      print(round(x$tdose, digits))
    } else if(!x$input$testOnly)
      cat("\nNone of significant models converged")
  }
  cat("\n")
}

summary.MCPMod <-
  function(object, digits = 3,...)
{
  oldClass(object) <- "summary.MCPMod"
  print(object, digits = digits)
}

print.summary.MCPMod <-
  function(x, digits = 3,...)
{
  cat("MCPMod\n\n")
  if(x$input$twoSide) side <- "two-sided"
  else side <- "one-sided"
  cat("Input parameters:","\n")
  if(attr(x$cVal, "Calc")){
    cat(" alpha =", paste(x$input$alpha," (", side,")", sep=""), "\n")
  }
  if(!x$input$testOnly){
    cat(" model selection:", paste(x$input$selModel[1]),"\n")
    if(is.element(x$input$selModel[1], c("aveAIC", "aveBIC"))){
      pW <- attr(x$fm, "pweights")
      cat(" prior model weights:\n ")
      print(round(pW/sum(pW), digits))
    }
    nr <- match(substr(x$input$doseEst,1,2), c("ME", "ED"))
    if(nr == 1){
      cat(" clinical relevance =",paste(x$input$clinRel),"\n")
    }
    dePar <- c("gamma", "p")[nr]
    cat(paste(" dose estimator: ", x$input$doseEst, " (",dePar, " = ", x$input$dePar, ")", sep=""), "\n")
  } # multiple contrast test information
  cat("\n","Optimal Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  cat("\n","Contrast Correlation:","\n", sep="")
  print(round(x$corMat, digits))
  cat("\n","Multiple Contrast Test:","\n",sep="")
  ord <- rev(order(x$tStat))
  if(x$input$pVal){
    print(data.frame(Tvalue = round(x$tStat, digits)[ord],
                     pValue = round(attr(x$tStat, "pVal"), digits)[ord]))
  }
  else {
    print(data.frame(Tvalue = round(x$tStat, digits)[ord]))
  }
  cat("\n","Critical value: ", round(x$cVal, digits),"\n",sep="")
  if(!x$signf) return(cat(""))  # No significant model
  if(all(is.na(x$model2))) {
    if(!x$input$testOnly) cat("\nNone of significant models converged","\n")
    return(cat(""))
  }
  else { # add IC values
    if(x$input$selModel[1] != "maxT"){
      if(x$input$selModel[1] == "AIC" | x$input$selModel[1]=="aveAIC")
        cat("\n",paste("AIC criterion:"),"\n",sep="")
      else cat("\n",paste("BIC criterion:"),"\n",sep="")
      print(round(attr(x$fm,"IC"), 2))
    } # model selection
    cat("\n","Selected for dose estimation:","\n", sep="")
    cat(" ", paste(x$model2),"\n\n", sep=" ")
    if(is.element(x$input$selModel[1], c("maxT", "AIC", "BIC"))){
      cat("Parameter estimates:","\n")
      cat(paste(x$model2), "model:\n")
      cof <- coef(x$fm[[1]])
      namcof <- names(cof)
      namcof <- gsub(" ", "", namcof)  # remove white spaces for GUI
      names(cof) <- gsub("doseM", "dose", namcof) # use more obvious names
      print(round(cof, digits))
    } else { # model averaging
      cat("Model weights:","\n", sep="")
      print(round(attr(x$fm, "weights"), digits))
      cat("\nParameter estimates:","\n")
      for(i in which(!is.na(attr(x$fm, c("IC"))))){
        nam <- names(x$fm)[i]
        cof <- coef(x$fm[[i]])
        namcof <- names(cof)
        namcof <- gsub(" ", "", namcof) # remove white spaces for GUI
        names(cof) <- gsub("doseM", "dose", namcof) # use more obvious names
        cat(paste(nam), "model:\n")
        print(round(cof, digits))
      }
    }
  } # information about dose estimate
  cat("\nDose estimate","\n")
  attr(x$tdose,"tdType") <- NULL # remove attr for output
  if(is.element(x$input$selModel, c("AIC", "BIC", "maxT"))) print(x$tdose)
  else {
    cat("Estimates for models\n")
    print(round(attr(x$tdose,"tdModels"), digits))
    attr(x$tdose,"tdModels") <- NULL
    cat("Model averaged dose estimate\n")

    print(round(x$tdose, digits))
  }
}

plot.MCPMod <-
  function(x, complData = FALSE, CI = FALSE, clinRel = FALSE, doseEst = FALSE, 
           gamma = NULL, models = "all", nrDoseGam = 1,
           colors = c("black","blue","black","gray","blue"),
           uGrad = NULL, ...)
{
  ## complData - logical indicating whether complete 
  ##             data set (T) or just group means 
  ##             should be plotted (F)
  ## CI - logical indicating whether confidence intervals 
  ##       should be plotted
  ## clinRel - logical indicating whether clinical relevance 
  ##      threshold should be plotted
  ## doseEst - logical whether dose estimate should be plotted
  ##       in case a vector was specified for dePar in the 
  ##       MCPMod function, nrDoseGam determines which is plotted
  ## gamma - value for the 1-2*gamma pointwise CI around
  ##         the predicted mean. if equal to NULL the value
  ##         determined in the MCPMod call is used. In case
  ##         a vector of gamma values nrDoseGam determines which
  ##         is used (see nrDoseGam)
  ## models - a character vector determining, which of the fitted
  ##          models should be plotted (only for model averaging)
  ## nrDoseGam - in case a vector is specified for gamma in 
  ##            the MCPMod function (and gamma in the plot.MCPMod
  ##            function is NULL), nrDoseGam determines which of
  ##            these values should be used for the conf. interval
  ##            and the dose estimate (if doseEst = T)
  ## colors - numeric of length 5 with number of colors for:
  ##          predictions, CI, data, clinical relevance threshold,
  ##          dose estimator
  ## ... - additional arguments for xyplot

  selModel <- x$input$selModel[1]
  ## model selection or model averaging?
  ms <- is.element(selModel, c("maxT", "AIC", "BIC")) 
  off <- x$input$off
  scal <- x$input$scal
  if(!is.null(x$input$resp)) resp <- x$input$resp # name of resp column
  else resp <- "resp"
  if(!is.null(x$input$dose)) dose <- x$input$dose # name of dose column
  else dose <- "dose"
  ind <- match(c(dose, resp), names(x$data))
  Data <- x$data[, ind]
  ## if no model is significant or no model converged plot raw data
  if(!x$signf | all(is.na(x$model2))){
    plot(Data[,1], Data[,2], xlab="Dose", ylab="Response")
    if(!x$signf) msg <- "No model significant"
    else msg <- "None of significant models converged"
    title(sub = msg)
    return(cat(msg,"\n"))
  } 
  
  if(is.null(gamma)){ # default for gamma
    if(x$input$doseEst != "ED"){
      gamma <- x$input$dePar[nrDoseGam] # use the same gamma as
    } else {                          # for dose estimate
      gamma <- 0.025 # if "ED" use reas. default
    }
  }
  if(models == "all"){ # select a subset of fitted models
    fm <- x$fm         # (only possible in case of model averaging)
    nam <- x$model2
  } else {
    if(ms) stop("models argument can only be used with model averaging")
    fm <- x$fm[models]
    nam <- models
  }
  lenDose <- x$input$lenDose # fit models and obtain CI
  doseSeq <- seq(min(Data[,1]), max(Data[,1]), length=lenDose)
  predList <- DoseInd <- list()
  z <- 0
  for(i in 1:length(fm)){
    if(is.list(fm[[i]])){
      z <- z + 1
      pred <- predict(fm[[i]], doseSeq, se.fit=CI, uGrad = uGrad)
      if(CI){
        LL <- pred$fit - qt(1-gamma, pred$df)*pred$se.fit
        UL <- pred$fit + qt(1-gamma, pred$df)*pred$se.fit
      }
      if(ms){
        tdose <- x$tdose[nrDoseGam] # plot just dose corresp. to
                                    # selected gamma value
      } else {
        tdose <- attr(x$tdose, "tdModels")
        modNum <- which(names(fm)[i]==dimnames(tdose)[[2]])
        tdose <- tdose[nrDoseGam, modNum]
      }
      DoseInd[[z]] <- ifelse(is.na(tdose), NA, which(doseSeq == tdose))
      if(CI) predList[[z]] <- c(pred$fit, LL, UL)
      else predList[[z]] <- pred
    }
  }
  plotVec <- do.call("c", predList)
  DoseInd <- do.call("c", DoseInd)
  if(CI) groups <- c("pred","LL","UL")
  else groups <- "pred"
  lg <- length(groups)
  plotMat <- data.frame(pred=plotVec, dose=rep(doseSeq, lg*z), 
                        nam=rep(nam, each=lg*lenDose),
                        group=rep(rep(groups, rep(lenDose, lg)), z))
  if(complData) { # plot all data, not just means
    dat <- Data
  } else {
    me <- tapply(Data[,2], Data[,1], mean)
    d <- as.numeric(names(me))    
    dat <- data.frame(d, me)
  }
  rg <- range(dat[,2], plotMat$pred)
  yl <- c(rg[1] - 0.1*diff(rg), rg[2] + 0.1*diff(rg))
  panDat <- list(data=dat, clinRel=x$input$clinRel, complData=complData,
                 cR=clinRel, doseEst=doseEst, lenDose=lenDose, DoseInd=DoseInd,
                 colors=colors, lg=lg) # data for panels
  ## information for key argument
  txt <- ifelse(complData, "Responses", "Group Means")
  pchar <- ifelse(complData, 1, 16)
  keyList <- list(points= list(col = colors[3], pch=pchar),
                  text = list(lab = txt),
                  lines = list(col = colors[1]),
                  text = list(lab = "Model Predictions"))
  if(CI){
    keyList <- c(keyList, list(lines = list(col = colors[2])), 
                 list(text = list(lab = paste(1-2*gamma, "Pointwise CI"))))
    col <- c(colors[2],colors[1],colors[2])
  } else col <- c(colors[1])
  if(doseEst){
    keyList <- c(keyList, list(points=list(col=colors[5], pch=18)),
                 list(text=list(lab="Estim. Dose")))
  }
  xyplot(pred~dose|nam, data=plotMat, xlab = "Dose", type="l",
         col=col, ylab = "Response",
         groups=plotMat$group, pD=panDat, ylim = yl, 
         panel=function(x, y, subscripts, groups,..., pD) {
           panel.superpose(x,y,subscripts,groups, ...)
           if(!pD$complData){
             panel.xyplot(pD$data[,1],pD$data[,2], pch=16,
                          col=pD$colors[3]) # plot data/means
           } else {
             panel.xyplot(pD$data[,1],pD$data[,2], col=pD$colors[3])
           }
           if(pD$cR) panel.abline(h=y[1]+pD$clinRel, lty=8,
                                  col=pD$colors[4]) # plot clinical
                                                    # relevance threshold
           if(pD$doseEst) {
             ind <- max(subscripts)/(pD$lenDose*pD$lg) # locate dose
                                                       # estimator in x
             panel.dotplot(x[pD$DoseInd[ind]], y[pD$DoseInd[ind]],
                           pch=18, col=pD$colors[5], col.line=0) # plot dose  
           }
         },
         strip = function(...) strip.default(..., style = 1), as.table=TRUE,
         key = keyList, ...)
}

#### example code
### 
### mods <- list(linear = NULL, linlog = NULL, emax = 0.3,
###          quadratic = -0.001, logistic = c(0.4, 0.09))
###
### mn <- c(0.1, 0.4, 0.55, 0.75, 0.9, 1)
### ds <- c(0, 0.05, 0.2, 0.4, 0.6, 1)
### DRdata <- genDFdata(mu = mn, sigma = 0.8, n = 20, doses = ds) 
### MM <- MCPMod(DRdata, mods, clinRel = 0.5*max(mn), selModel="aveAIC")
### MM
### summary(MM)
### plot(MM)
###
### user defined model
### emx2 <- function(dose, a, b, d){
###   sigEmax(dose, a, b, d, 1)
### }
### emx2Grad <- function(dose, a, b, d) cbind(1, dose/(dose+d), -b*dose/(dose+d)^2)
### models <- list(emx2=c(0,1,0.1), linear = NULL)
### dats <- genDFdata(mu = mn, doses = ds, sigma=0.8, n=50)
### start <- list(emx2=c(a=0.2, b=0.6, d=0.2))
### MM1 <- MCPMod(dats, models, clinRel = 1, scal = 1.2, selModel="AIC", start = start,
###          uGrad = emx2Grad)

## function to generate DF data
genDFdata <-
  function(model, argsMod, doses, n, sigma, mu = NULL)
{
  ## generates data.frame with resp and dose columns
  ## corresponding to specified design, mean (either
  ## determined by model or mu) and sigma
  ##
  ## model - string with model function name (first argument
  ##         must be dose, must be vectorized)             
  ## argsMod - a named list with the arguments for the model function
  ##
  ## doses - set of doses to be used in the trial
  ## n - sample size (scalar is repeated for each dose)
  ## sigma - error std. deviation
  ## mu - vector of mean values can alternatively specified

  nD <- length(doses)
  dose <- sort(doses)
  if (length(n) == 1) n <- rep(n, nD)
  dose <- rep(dose,  n)
  if(!missing(model)){
    args <- c(list(dose), argsMod)
    mu <- do.call(model, args)
  } else if(!is.null(mu)){
      if(length(doses) != length(mu)){
        stop("'mu' needs to be of the same length as doses")
       }
       mu <- rep(mu,  n)
    } else {
    stop("either 'model' or 'mu' needs to be specified")
  }
  data.frame(dose = dose, 
             resp = mu + rnorm(sum(n), sd = sigma))
}
### examples
### genDFdata("emax", c(e0 = 0.2, eMax = 1, ed50 = 0.05), c(0,0.05,0.2,0.6,1), 20, 1)
### genDFdata(mu = 1:5, doses = 0:4, n = c(20, 20, 10, 5, 1), sigma = 1)
