# file mdscore/R/mdscore.R
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

mdscore <- function(model = model, X1 = X1, phi = NULL){

  arg.phi <- phi  
  if(!is.element("glm", class(model))) 
    stop(paste("\n Model should be a 'glm' object", "!!!\n"))
                                                                       
  fam <- summary(model)$family$family 
  if(is.null(arg.phi) & !(fam %in% c("gaussian", "Gamma", "inverse.gaussian", "poisson", "binomial"))) 
    stop(paste("\n When the parameter phi is unknown, the model family must be 'gaussian', 'Gamma' or 'inverse.gaussian'","!!!\n",
               "\n When the parameter phi is fixed, the model family must be 'poisson' or 'binomial'", "!!!\n"
               )
         )
  
  X1 <- as.matrix(X1)
  X2 <- model.matrix(model)
  X  <- cbind(X1,X2)
  n  <- dim(X)[1]
  p  <- dim(X)[2]
  q  <- dim(X1)[2]
  
  #----------------------------------------------------------------------
  # Mean and variance functions
  #----------------------------------------------------------------------
  if(model$family[[2]] %in% c("log", "cloglog", "logit")){
    
    mu <- switch(model$family[[2]],
                 log     = quote(exp(eta)),
                 cloglog = quote(1 - exp(-exp(eta))),
                 logit   = quote(exp(eta)/(1 + exp(eta))))
    
  }else mu <- as.list(model$family$linkinv)[[2]]
       
  V <- if(model$family[[1]] == "gaussian") quote(1) else
       as.list(model$family$variance)[[2]]
    
  #----------------------------------------------------------------------
  # Mean and variance derivatives
  #----------------------------------------------------------------------
  dmu  <- D(mu,"eta")
  d2mu <- D(dmu,"eta")
  dV   <- D(V,"mu")
  d2V  <- D(dV,"mu")
  
  #----------------------------------------------------------------------
  # Parameter estimates
  #----------------------------------------------------------------------
  mu.est   <- model$fitted.values
  eta.est  <- model$family$linkfun(mu.est)
  V.est    <- eval(V,list(mu = mu.est))
  dmu.est  <- as.numeric(eval(dmu,list(eta = eta.est)))
  d2mu.est <- as.numeric(eval(d2mu,list(eta = eta.est)))
  dV.est   <- as.numeric(eval(dV,list(mu = mu.est)))
  d2V.est  <- as.numeric(eval(d2V,list(mu = mu.est)))  
  if(is.null(arg.phi)) phi <- estphi(model)
  
  #----------------------------------------------------------------------
  # Vectors and matrices
  #----------------------------------------------------------------------  
  f <- 1/V.est*dmu.est*d2mu.est
  g <- f - 1/V.est^2*dV.est*(dmu.est)^3
  b <- 1/V.est^3*dmu.est^4*(dV.est^2 + V.est*d2V.est)
  h <- 1/V.est^2*(dV.est*dmu.est^2*d2mu.est + d2V.est*dmu.est^4)
  if(length(f) == 1) f <- rep(f,n)
  if(length(g) == 1) g <- rep(g,n)
  if(length(b) == 1) b <- rep(b,n)
  if(length(h) == 1) h <- rep(h,n)
  
  w <- as.numeric(model$weights) # w = dmu.est^2/V.est  
  IX2WX2 <- invmat(crossprod(X2, w*X2)) 
  X2WX1  <-  crossprod(X2, w*X1)  
  C      <- crossprod(IX2WX2, X2WX1) 
  R      <- X1 - crossprod(t(X2), C) 
  IRWR   <- invmat(crossprod(R, w*R)) 
  IXWX   <- rbind(cbind(IRWR, -IRWR%*%t(C)), 
                  cbind(-C%*%IRWR, IX2WX2 + C%*%IRWR%*%t(C))) 
  Z      <- crossprod(t(X), tcrossprod(IXWX, X)) 
  Z2     <- crossprod(t(X2), tcrossprod(IX2WX2, X2)) 
  Z_Z2   <- crossprod(t(R), tcrossprod(IRWR, R))
  Z2d       <- diag(diag(Z2))
  Zd        <- diag(diag(Z))
  ZminusZ2d <- diag(diag(Z_Z2))
  
  #----------------------------------------------------------------------
  # Correction terms
  #----------------------------------------------------------------------  
  i1 <- rep(1,n)  
  A1.beta <- phi^-1*(+3*t(i1)%*%(f*Z2d)%*%Z_Z2%*%(Z2d*f)%*%i1
                     +6*t(i1)%*%(f*Z2d)%*%Z2%*%(ZminusZ2d*(f-g))%*%i1
                     -6*t(i1)%*%(f*(Z2^2)*((Z_Z2))*(2*g-f))%*%i1
                     -6*t(i1)%*%(h*ZminusZ2d)%*%Z2d%*%i1)
  A2.beta <- phi^-1*(-3*t(i1)%*%((f-g)*ZminusZ2d)%*%Z2%*%(ZminusZ2d*(f-g))%*%(i1)
                     -6*t(i1)%*%(f*Z2d)%*%(Z_Z2)%*%(ZminusZ2d*(f-g))%*%i1
                     -6*t(i1)%*%(((f-g)*((Z_Z2)^2))*(Z2*(f-g)))%*%i1
                     +3*t(i1)%*%(b*ZminusZ2d^2)%*%i1)
  A3.beta <- phi^-1*(+3*t(i1)%*%((f-g)*ZminusZ2d)%*%(Z_Z2)%*%(ZminusZ2d*(f-g))%*%i1
                     +2*t(i1)%*%(((f-g)*((Z_Z2)^3))*(f-g))%*%i1)
  
  if(is.null(arg.phi)){
    
    if(model$family$family == "Gamma"){
      
      A1.beta.phi <- -6*q*(1 + phi^2*psigamma(phi,deriv=2) - (2 - p + q)*(1 - phi*trigamma(phi)))/
                      (n*phi*(1-phi*trigamma(phi))^2)
      
      A2.beta.phi <- (3*q*(q + 2))/(n*phi*(1 - phi*trigamma(phi)))
      
    } else{
      
      if(model$family$family == "gaussian" | model$family$family == "inverse.gaussian"){
        
        A1.beta.phi <- 12*q*(p - q)/n
        A2.beta.phi <- -6*q*(q + 2)/n
        
      } else{
        
        A1.beta.phi <- 0
        A2.beta.phi <- 0
        
      }
      
    }
    
  } else{
    
    A1.beta.phi <- 0
    A2.beta.phi <- 0
    
  }
  
  A1 <- A1.beta + A1.beta.phi
  A2 <- A2.beta + A2.beta.phi
  A3 <- A3.beta
  
  a <- as.numeric(A3/(12*q*(q + 2)*(q + 4)))
  b <- as.numeric((A2 - 2*A3)/(12*q*(q + 2)))
  c <- as.numeric((A1 - A2 + A3)/(12*q))
  
  #----------------------------------------------------------------------
  # Pearson residuals
  #----------------------------------------------------------------------  
  rp <- resid(model,type="pearson")*sqrt(phi) #also rp <- sqrt(phi)*(model$y - mu.est)/sqrt(V.est)
  
  #----------------------------------------------------------------------
  # Score and corrected score statistic
  #----------------------------------------------------------------------  
  Sr     <- as.numeric(t(rp)%*%(sqrt(w)*X1)%*%IRWR%*%t(X1*sqrt(w))%*%rp)
  Sr_cor <- as.numeric(Sr*(1-(c+b*Sr+a*Sr^2)))  
  if(Sr_cor <= 0) stop("Correction yields a non positive value") 
  
#  return(list(Sr=Sr,Sr_cor=Sr_cor,coef=c(a,b,c),n=n,df=q, phi=phi))

  structure(
         list(
            Sr = Sr,
            Sr_cor = Sr_cor,
            coef = c(a,b,c),
            n = n,
            df = q,
            phi = phi
         ),
         class="mdscore"
  )
}

# Function estphi: Estimates the phi parameter 
estphi <- function(m){ 
  
  n <- nrow(model.matrix(m))
  k <- ncol(model.matrix(m))
  switch(summary(m)$family$family, 
         gaussian         = 1/summary(m)$dispersion*n/summary(m)$df.residual, 
         Gamma            = 1/gamma.dispersion(m), 
         inverse.gaussian = 1/summary(m)$dispersion*n/summary(m)$df.residual,
         poisson          = summary(m)$dispersion,
         binomial         = summary(m)$dispersion
  )
}

#----------------------------------------------------------------------
# Routine for matrix inversion
#----------------------------------------------------------------------

invmat <- function(mat){
  
  if(nrow(mat) <= 1) return(1/mat) 
  else return(chol2inv(chol(mat))) 
  
}

summary.mdscore <- function(object, ...){
  if (!inherits(object, "mdscore")) 
    warning("calling summary.mdscore(<fake-mdscore-object>) ...")
  
  stats <- round(c(object$Sr, object$Sr_cor), 2)
  pvals <- c(pchisq(object$Sr, df=object$df, lower.tail=FALSE),
             pchisq(object$Sr_cor, df=object$df, lower.tail=FALSE))
  pvals <- ifelse(pvals>=0.0001, round(pvals, 4), "< 0.0001")           
  tab   <- data.frame(object$df, stats, pvals)
  
  rownames(tab) <- c("Score", "Modified score")
  names(tab)    <- c("Df", " Value", " P-value")
  cat("\n\n")
  return(tab)
}