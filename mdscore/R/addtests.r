# file mdscore/R/addtests.R
# by  Damiao N. da Silva and Antonio Hermes M. da Silva Jr.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

wald.gamma <- function(model = model, terms){
  
  X   <- model.matrix(model)
  X1  <- as.matrix(X[,terms])
  X2  <- as.matrix(X[,-terms])
  w   <- as.numeric(model$weights)
  pp  <- as.numeric(gamma.shape(model)$alpha) 
  
  IX2WX2 <- invmat(crossprod(X2, w*X2))  
  X2WX1  <- crossprod(X2, w*X1)  
  C      <- crossprod(IX2WX2, X2WX1) 
  R      <- X1 - crossprod(t(X2), C) 
  RWR    <- crossprod(R, w*R) 
  b1     <- coef(model)[terms]
  
  est   <- pp*t(b1)%*%RWR%*%b1
  pvalue <- pchisq(est, df=ncol(X1), lower.tail=FALSE)
  list(est=est, pvalue=pvalue)
}

lrt.gamma <- function(y, fit1, fit2){
  
  phi1 <- gamma.shape(fit1)$alpha
  phi2 <- gamma.shape(fit2)$alpha
  
  S1 <- 2*sum(dgamma(y, shape=phi1, scale=fit1$fitted/phi1, log = TRUE))
  S2 <- 2*sum(dgamma(y, shape=phi2, scale=fit2$fitted/phi2, log = TRUE))
  est <- S2 - S1 
  
  g1 <- fit2$df.residual
  g2 <- fit1$df.residual
  
  if(est < 0){
    est <- abs(est)
    warning(paste("\n 'fit1' must be under null hypothesis","!!!\n"))
    g2 <- fit1$df.residual
    g1 <- fit2$df.residual
  }
  
  pvalue <- 1 - pchisq(est, g2 - g1)
  list(est=est, pvalue=pvalue)
}


wald.norm <- function(model = model, terms){
  
  X   <- model.matrix(model)
  X1  <- as.matrix(X[,terms])
  X2  <- as.matrix(X[,-terms])
  w   <- as.numeric(model$weights)
  n   <- nrow(X) 
  pp  <- as.numeric(1/summary(model)$dispersion)*n/summary(model)$df.residual
  
  IX2WX2 <- invmat(crossprod(X2, w*X2))  
  X2WX1  <- crossprod(X2, w*X1)  
  C      <- crossprod(IX2WX2, X2WX1) 
  R      <- X1 - crossprod(t(X2), C) 
  RWR    <- crossprod(R, w*R) 
  b1     <- coef(model)[terms]
  
  est    <- pp*t(b1)%*%RWR%*%b1
  pvalue <- pchisq(est, df=ncol(X1), lower.tail=FALSE)
  list(est=est, pvalue=pvalue)  
  
}

lrt.norm <- function(y,fit1,fit2){
  
  n    <- length(fit1$fitted)
  phi1 <- n/fit1$deviance
  phi2 <- n/fit2$deviance
  
  S1  <- 2*sum(dnorm(y, mean=fit1$fitted, sd=1/sqrt(phi1), log = TRUE))
  S2  <- 2*sum(dnorm(y, mean=fit2$fitted, sd=1/sqrt(phi2), log = TRUE))
  est <- S2 - S1  
  
  g1 <- fit2$df.residual
  g2 <- fit1$df.residual
  
  if(est < 0){
    est <- abs(est)
    warning(paste("\n 'fit1' must be under null hypothesis","!!!\n"))
    g2 <- fit1$df.residual
    g1 <- fit2$df.residual
  }
  
  pvalue <- 1 - pchisq(est, g2 - g1)
  list(est=est, pvalue=pvalue)
}

wald.invg <- function(model = model, terms){
  
  X   <- model.matrix(model)
  X1  <- as.matrix(X[,terms])
  X2  <- as.matrix(X[,-terms])
  w   <- as.numeric(model$weights)
  n   <- nrow(X) 
  pp  <- as.numeric(1/summary(model)$dispersion)*n/summary(model)$df.residual
  
  IX2WX2 <- invmat(crossprod(X2, w*X2))  
  X2WX1  <- crossprod(X2, w*X1)  
  C      <- crossprod(IX2WX2, X2WX1) 
  R      <- X1 - crossprod(t(X2), C) 
  RWR    <- crossprod(R, w*R) 
  b1     <- coef(model)[terms]
  
  est    <- pp*t(b1)%*%RWR%*%b1
  pvalue <- pchisq(est, df=ncol(X1), lower.tail=FALSE)
  return(list(est=est, pvalue=pvalue))
}

dinvg <- function(x, mu, phi, log=FALSE){
  
  if(any(x<=0))  stop("Values of x must be positive")
  if(any(mu<=0)) stop("Values of parameter mu must be positive")
  if(phi<=0)     stop("Value of parameter phi must be positive")
  
  dens <- -0.5*phi*(x - mu)^2/(x*mu^2) + 0.5*log(phi/(2*pi*x^3))
  if(log) return(dens) else return(exp(dens))
  
}

lrt.invg <- function(y, fit1, fit2){
  
  n    <- length(y)
  phi1 <- n/fit1$deviance
  phi2 <- n/fit2$deviance
  
  S1  <- 2*sum(dinvg(y, mu=fit1$fitted, phi=phi1, log = TRUE))
  S2  <- 2*sum(dinvg(y, mu=fit2$fitted, phi=phi2, log = TRUE))
  est <- S2 - S1 
  
  g1 <- fit2$df.residual
  g2 <- fit1$df.residual

  if(est < 0){
    est <- abs(est)
    warning(paste("\n 'fit1' must be under null hypothesis","!!!\n"))
    g2 <- fit1$df.residual
    g1 <- fit2$df.residual
  }
    
  pvalue <- 1 - pchisq(est, g2 - g1)
  list(est=est, pvalue=pvalue)  
}

wald.test <- function(model = model, terms){
  if(class(model)[1]!="glm")
    stop(paste("\n 'model' is not an object from 'glm' class", "!!!\n"))
  
  fam <- summary(model)$family$family 
  if(!(fam %in% c("gaussian", "Gamma", "inverse.gaussian"))) 
    stop(paste("\n When the parameter phi is unknown, the model family must be 'gaussian', 'Gamma' or 'inverse.gaussian'","!!!\n"))
  
  if(fam == "gaussian"){
    rslt = wald.norm(model = model, terms = terms)
  }else{
    if(fam == "Gamma"){
      rslt = wald.gamma(model = model, terms)
      }else{
        rslt = wald.invg(model = model, terms)
      }
    }
  structure(
    list(
      W = rslt$est,
      pvalue = rslt$pvalue
      ),
    class="wald.test"
    )   
}

lr.test <- function(fit1, fit2){
  if(class(fit1)[1]!="glm")
    stop(paste("\n 'fit1' is not an object from 'glm' class", "!!!\n"))
  
  if(class(fit2)[1]!="glm")
    stop(paste("\n 'fit2' is not an object from 'glm' class", "!!!\n"))
  
  if(!(all(fit1$y %in% fit2$y)))
    stop(paste("\n 'fit1' and 'fit2' must have the same response variable","!!!\n"))
    
  if(!(all(fit1$family %in% fit2$family)))
    stop(paste("\n 'fit1' and 'fit2' must have the same family and link function","!!!\n"))
  
  fam <- summary(fit1)$family$family 
  if(!(fam %in% c("gaussian", "Gamma", "inverse.gaussian"))) 
    stop(paste("\n When the parameter phi is unknown, the model family must be 'gaussian', 'Gamma' or 'inverse.gaussian'","!!!\n"))
  
  y = fit1$y

  if(fam == "gaussian"){
    rslt = lrt.norm(y, fit1, fit2)
  } else{
    if(fam == "Gamma"){
      rslt = lrt.gamma(y, fit1, fit2)
      } else{
        rslt = lrt.invg(y, fit1, fit2)
    }
  }
  
  structure(
    list(
      LR = rslt$est,
      pvalue = rslt$pvalue
    ),
    class="lrt.test"
  )   
}