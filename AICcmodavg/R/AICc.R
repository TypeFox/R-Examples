##compute AIC, AICc, QAIC, QAICc
##generic
AICc <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
  UseMethod("AICc", mod)
}

AICc.default <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
  stop("\nFunction not yet defined for this object class\n")
}



##aov objects
AICc.aov <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##betareg objects
AICc.betareg <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##clm objects
AICc.clm <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(fitted(mod))} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc       
  }


##clmm objects
AICc.clmm <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(fitted(mod))} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##coxme objects
AICc.coxme <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
  
  if(identical(nobs, NULL)) {n <- length(mod$linear.predictor)} else {n <- nobs}
  LL <- extractLL(mod)[1]
  K <- attr(extractLL(mod), "df")  #extract correct number of parameters included in model
  if(second.ord==TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
  if(return.K==TRUE) AICc[1] <- K #attributes the first element of AICc to K
  AICc
}



##coxph objects
AICc.coxph <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
  if(identical(nobs, NULL)) {n <- length(residuals(mod))} else {n <- nobs}
  LL <- extractLL(mod)[1]
  K <- attr(extractLL(mod), "df")  #extract correct number of parameters included in model
  if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
  if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
  AICc
}



##fitdist (from fitdistrplus)
AICc.fitdist <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- mod$n} else {n <- nobs}
    LL <- logLik(mod)
    K <- length(mod$estimate)
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##fitdistr (from MASS)
AICc.fitdistr <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- mod$n} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##glm and lm objects
AICc.glm <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1, ...){
    
    if(is.null(nobs)) {
      n <- length(mod$fitted)
    } else {n <- nobs}
    
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
      
    if(c.hat == 1) {
      if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))} else{AICc <- -2*LL+2*K}
    }
    if(c.hat > 1 && c.hat <= 4) {
      K <- K+1
      if(second.ord==TRUE) {
        AICc <- (-2*LL/c.hat)+2*K*(n/(n-K-1))
        ##adjust parameter count to include estimation of dispersion parameter
      } else{
        AICc <- (-2*LL/c.hat)+2*K}
    }

    if(c.hat > 4) stop("High overdispersion and model fit is questionable\n")
    if(c.hat < 1) stop("You should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

    ##check if negative binomial and add 1 to K for estimation of theta if glm( ) was used
    if(!is.na(charmatch(x="Negative Binomial", table=family(mod)$family))) {
      if(!identical(class(mod)[1], "negbin")) { #if not negbin, add + 1 because k of negbin was estimated glm.convert( ) screws up logLik
        K <- K+1
        if(second.ord == TRUE) {
          AICc <- -2*LL+2*K*(n/(n-K-1))
        } else {
          AICc <- -2*LL+2*K
        }
      }
      if(c.hat != 1) stop("You should not use the c.hat argument with the negative binomial")
    }
    ##add 1 for theta parameter in negative binomial fit with glm( )

    ##check if gamma and add 1 to K for estimation of shape parameter if glm( ) was used
    if(identical(family(mod)$family, "Gamma") && c.hat > 1) stop("You should not use the c.hat argument with the gamma")
      
    ##an extra condition must be added to avoid adding a parameter for theta with negative binomial when glm.nb( ) is fit which estimates the correct number of parameters
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##gls objects
AICc.gls <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n<-length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##gnls objects
AICc.gnls <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(fitted(mod))} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##hurdle objects
AICc.hurdle <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL + 2*K*(n/(n-K-1))}  else{AICc <- -2*LL + 2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }


##lavaan
AICc.lavaan <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- mod@Data@nobs[[1]]} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }


##lm objects
AICc.lm <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##lme objects
AICc.lme <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- nrow(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }


##lmekin objects
AICc.lmekin <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(mod$residuals)} else {n <- nobs}
    LL <- extractLL(mod)[1]
    K <- attr(extractLL(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    return(AICc)
  }


##maxlike objects
AICc.maxlikeFit <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1, ...) {
  
  LL <- extractLL(mod)[1]
  K <- attr(extractLL(mod), "df")

  if(is.null(nobs)) {
    n <- nrow(mod$points.retained)
  } else {n <- nobs}
  
  if(second.ord == TRUE) {AICc <- -2 * LL + 2 * K * (n/(n - K - 1))} else {AICc <- -2*LL + 2*K}
  
  if(c.hat != 1) stop("\nThis function does not support overdispersion in \'maxlikeFit\' models\n")

  if(identical(return.K, TRUE)) {
    return(K)
  } else {return(AICc)}
}


##mer object
AICc.mer <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(is.null(nobs)) {
      n <- mod@dims["n"]
    } else {n <- nobs}
      
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
    
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))} else{AICc <- -2*LL+2*K}
    
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }


##merMod objects
AICc.merMod <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
 
    if(is.null(nobs)) {
      n <- mod@devcomp$dims["n"]
      names(n) <- NULL
    } else {n <- nobs}
    
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
      
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))} else{AICc <- -2*LL+2*K}
                    
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }


##mult objects
AICc.multinom <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1, ...){

    if(identical(nobs, NULL)) {n<-length(mod$fitted)/length(mod$lev)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
      
    if(c.hat == 1) {
      if(second.ord==TRUE) {
        AICc <- -2*LL+2*K*(n/(n-K-1))
      } else{
        AICc <- -2*LL+2*K
      }
    }
    if(c.hat > 1 && c.hat <= 4) {
      K <- K+1
      if(second.ord == TRUE) {
        AICc <- (-2*LL/c.hat)+2*K*(n/(n-K-1)) #adjust parameter count to include estimation of dispersion parameter
      } else{
        AICc <- (-2*LL/c.hat)+2*K
      }
    }
    if(c.hat > 4) stop("High overdispersion and model fit is questionable")
    
    if(return.K==TRUE) AICc[1]<-K #attributes the first element of AICc to K
    AICc
  }


##nlme objects
AICc.nlme <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- nrow(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }


##nls objects
AICc.nls <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(fitted(mod))} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }


##polr objects
AICc.polr <-
function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

  if(identical(nobs, NULL)) {n<-length(mod$fitted)} else {n <- nobs}
  LL <- logLik(mod)[1]
  K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
      
  if(second.ord == TRUE) {
    AICc <- -2*LL+2*K*(n/(n-K-1))
  } else{
    AICc <- -2*LL+2*K
  }
                                 
  if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
  AICc
}



##rlm objects
##only valid for M-estimation (Huber M-estimator)
##modified from Tharmaratnam and Claeskens 2013 (equation 8)
##AICc.rlm <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...)
##{

##  if(second.ord == TRUE) stop("\nOnly 'second.ord = FALSE' is supported for 'rlm' models\n")

##  ##extract design matrix
##  X <- model.matrix(mod)
  
##  ##extract scale
##  scale.m <- mod$s

##  ##extract threshold value
##  cval <- mod$k2
    
##  ##extract residuals
##  res <- residuals(mod)
##  res.scaled <- res/scale.m
##  n <- length(res)
  
##  ##partial derivatives based on Huber's loss function
##  dPsi <- ifelse(abs(res.scaled) <= cval, 2, 0)
##  Psi <- (ifelse(abs(res.scaled) <= cval, 2*res.scaled, 2*cval*sign(res.scaled)))^2
    
##  J <- (t(X) %*% diag(as.vector(dPsi)) %*% X * (1/(scale.m^2)))/n
##  inv.J <- solve(J)
  
##  ##variance
##  K.var <- (t(X) %*% diag(as.vector(Psi)) %*% X * (1/(scale.m^2)))/n
##  AIC <- 2*n*(log(scale.m)) + 2 * sum(diag(inv.J %*%(K.var)))
  
##  if(return.K) {AIC <- 2 * sum(diag(inv.J %*%(K.var)))}
##  return(AIC)
##}

##the estimator below extracts the estimates obtained from M- or MM-estimator
##and plugs them in the normal likelihood function
AICc.rlm <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##survreg objects
AICc.survreg <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- nrow(mod$y)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##unmarkedFit objects
##create function to extract AICc from 'unmarkedFit'
AICc.unmarkedFit <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1, ...) {
  
  LL <- extractLL(mod)[1]
  K <- attr(extractLL(mod), "df")
  
  if(is.null(nobs)) {
    n <- dim(mod@data@y)[1]
  } else {n <- nobs}
  
  if(c.hat == 1) {
    if(second.ord == TRUE) {AICc <- -2 * LL + 2 * K * (n/(n - K - 1))} else {AICc <- -2*LL + 2*K}
  }
  if(c.hat > 1 && c.hat <= 4) {
    ##adjust parameter count to include estimation of dispersion parameter
    K <- K + 1
    if(second.ord == TRUE) {
      AICc <- (-2 * LL/c.hat) + 2 * K * (n/(n - K - 1))
    } else {
      AICc <- (-2 * LL/c.hat) + 2*K}
  }

  if(c.hat > 4) stop("\nHigh overdispersion and model fit is questionable\n")
  if(c.hat < 1) stop("\nYou should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

  if(identical(return.K, TRUE)) {
    return(K)
  } else {return(AICc)}
}


##vglm objects
AICc.vglm <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1, ...){
    
    if(is.null(nobs)) {
      n <- nrow(mod@fitted.values)
    } else {n <- nobs}
    
    LL <- extractLL(mod)[1]

    ##extract number of estimated parameters
    K <- attr(extractLL(mod), "df")
    
    if(c.hat !=1) {
      fam.name <- mod@family@vfamily
      if(fam.name != "poissonff" && fam.name != "binomialff") stop("\nOverdispersion correction only appropriate for Poisson or binomial models\n")
    }
    if(c.hat == 1) {
      if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))} else{AICc <- -2*LL+2*K}
    }
    if(c.hat > 1 && c.hat <= 4) {
      K <- K + 1
      if(second.ord==TRUE) {
        AICc <- (-2*LL/c.hat) + 2*K*(n/(n-K-1))
        ##adjust parameter count to include estimation of dispersion parameter
      } else{
        AICc <- (-2*LL/c.hat) + 2*K}
    }

    if(c.hat > 4) stop("High overdispersion and model fit is questionable\n")
    if(c.hat < 1) stop("You should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }



##zeroinfl objects
AICc.zeroinfl <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL + 2*K*(n/(n-K-1))}  else{AICc <- -2*LL + 2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  }

