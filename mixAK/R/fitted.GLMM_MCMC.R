##
##  PURPOSE:  Function to calculate fitted profiles from GLMM_MCMC object
##            based on posterior means of regression coefficients
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   24/02/2010 (as a stand alone function)
##             26/11/2010:  added to the mixAK package
##             27/12/2010:  argument statistic added
##             17/09/2012:  warning for code under overall=TRUE added
##             22/09/2013:  Gaussian quadrature provided by package fastGHQuad and no more
##                          by logLik from glmer
##
##  FUNCTIONS: fitted.GLMM_MCMC
##
## ==========================================================================

## *************************************************************
## fitted.GLMM_MCMC
## *************************************************************
##
fitted.GLMM_MCMC <- function(object, x, z, statistic=c("median", "mean", "Q1", "Q3", "2.5%", "97.5%"), overall=FALSE,
                             glmer=TRUE, nAGQ=100, ...)
{
  if (overall){
    warning("with overall=FALSE this function might be imprecise.")
    cat("The fitted longitudinal profiles for response variables are calculated\n")
    cat("as if the random effects are normally distributed with the mean given\n")
    cat("by the posterior summary statistic of the random effects overall mean\n")
    cat("and the covariance matrix given by the posterior summary statistic\n")
    cat("of the random effects overall covariance matrix.\n")
    cat("That is, the original normal mixture is approximated\n")
    cat("by a one-component normal distribution.\n")
  }  

  QRule <- fastGHQuad::gaussHermiteData(nAGQ)       ### needed by fastGHQuad::aghQuad function used in case of a binomial response
  
  statistic <- match.arg(statistic)
  
  if (statistic == "median") Stat <- "Median"
  else if (statistic == "mean") Stat <- "Mean"
       else if (statistic == "Q1") Stat <- "1st Qu."
            else if (statistic == "Q3") Stat <- "3rd Qu."
                 else if (statistic == "2.5%") Stat <- "2.5%"
                      else if (statistic == "97.5%") Stat <- "97.5%"
  
  RR <- sum(object$R)
  
  if (object$prior.b$priorK != "fixed") stop("implemented only for models with a fixed number of mixture components")
  K_b <- object$K_b[1]
  
  if (sum(object$p)){
    if (missing(x)) stop("x must be given")    
    if (RR == 1 & !is.list(x)) x <- list(x)
    if (!is.list(x)) stop("x must be a list")
    if (length(x) != RR) stop("x must be a list of length", RR)
  }
  
  if (sum(object$q)){
    if (missing(z)) stop("z must be given")
    if (RR == 1 & !is.list(z)) z <- list(z)
    if (!is.list(z)) stop("z must be a list")
    if (length(z) != RR) stop("z must be a list of length", RR)
  }    
  
  #### Posterior summary statistic of means of random effects in each component
  #### and possibly also for covariance matrices
  if (object$dimb){
    if (overall){
      K_b <- 1
      if (object$dimb == 1) mu_b <- matrix(object$summ.b.Mean[Stat], nrow=1)
      else                  mu_b <- matrix(object$summ.b.Mean[Stat,], nrow=1)

      if (glmer){
        if (object$dimb == 1) D_b <- matrix(object$summ.b.SDCorr[Stat]^2, nrow=1, ncol=1)
        else{
          D_b <- matrix(NA, nrow=object$dimb, ncol=object$dimb)
          for (j in 1:object$dimb){
            D_b[j, j] <- object$summ.b.SDCorr[Stat, paste("b.SD.", j, sep="")]^2
            if (j == object$dimb) break
            for (i in (j+1):object$dimb) D_b[i, j] <- D_b[j, i] <- object$summ.b.SDCorr[Stat, paste("b.Corr.", i, ".", j, sep="")] * object$summ.b.SDCorr[Stat, paste("b.SD.", i, sep="")] * object$summ.b.SDCorr[Stat, paste("b.SD.", j, sep="")]
          }  
        }
        D_b <- list(D_b)
      }  
      
    }else{
      mu_b <- object$poster.mean.mu_b * matrix(object$scale.b$scale, nrow=K_b, ncol=object$dimb, byrow=TRUE) + matrix(object$scale.b$shift, nrow=K_b, ncol=object$dimb, byrow=TRUE)

      if (glmer){
        D_b <- list()
        for (k in 1:K_b){
          if (object$dimb == 1) D_b[[k]] <- object$scale.b$scale^2 * object$poster.mean.Sigma_b[[k]]
          else                  D_b[[k]] <- diag(object$scale.b$scale) %*% object$poster.mean.Sigma_b[[k]] %*% diag(object$scale.b$scale)
        }  
      }        
    }  
  }

  #### Posterior summary statistic of fixed effects
  if (object$lalpha){
    if (object$lalpha == 1) alpha <- object$summ.alpha[Stat]
    else{
      alpha <- as.numeric(object$summ.alpha[Stat,])
      names(alpha) <- colnames(object$summ.alpha)
    }  
  }  
  
  #### Loop over response types
  qri <- object$q + object$random.intercept
  pfi <- object$p + object$fixed.intercept

  fit <- mfit <- list()
  for (r in 1:RR){        ### loop (rr)
    
    ### there are random effects
    if (qri[r]){
      
      if (r == 1){
        Eb <- if (K_b > 1) matrix(mu_b[, 1:qri[r]], nrow=K_b) else matrix(mu_b[, 1:qri[r]], nrow=1)

        if (glmer){
          Vb <- list()
          if (qri[r] > 1) for (k in 1:K_b) Vb[[k]] <- D_b[[k]][1:qri[r], 1:qri[r]]
          else            for (k in 1:K_b) Vb[[k]] <- matrix(D_b[[k]][1:qri[r], 1:qri[r]], nrow=1)
        }        
      }else{  
        Eb <- if (K_b > 1) matrix(mu_b[, (sum(qri[1:(r-1)])+1):sum(qri[1:r])], nrow=K_b) else matrix(mu_b[, (sum(qri[1:(r-1)])+1):sum(qri[1:r])], nrow=1)

        if (glmer){
          Vb <- list()
          if (qri[r] > 1) for (k in 1:K_b) Vb[[k]] <- D_b[[k]][(sum(qri[1:(r-1)])+1):sum(qri[1:r]), (sum(qri[1:(r-1)])+1):sum(qri[1:r])]
          else            for (k in 1:K_b) Vb[[k]] <- matrix(D_b[[k]][(sum(qri[1:(r-1)])+1):sum(qri[1:r]), (sum(qri[1:(r-1)])+1):sum(qri[1:r])], nrow=1)
        }  
      }  

      if (object$q[r]){
        if (object$q[r] == 1){
          if (is.matrix(z[[r]])) if (ncol(z[[r]]) != 1) stop("z[[", r, "]] must have 1 column", sep="") 
          z[[r]] <- matrix(z[[r]], ncol=1)
        }
        
        if (!is.matrix(z[[r]])) stop("z[[", r, "]] must be a matrix", sep="")
        if (ncol(z[[r]]) != object$q[r]) stop("z[[", r, "]] must have ", object$q[r], " columns", sep="")
      }  

      if (object$random.intercept[r]){
        if (object$q[r]){
          fit[[r]] <- Eb[1, 1] + z[[r]] %*% Eb[1, 2:qri[r]]
          if (K_b > 1) for (k in 2:K_b) fit[[r]] <- cbind(fit[[r]], Eb[k, 1] + z[[r]] %*% Eb[k, 2:qri[r]])
        }else{
          fit[[r]] <- matrix(Eb, nrow=1, ncol=K_b)
        }          
      }else{
        fit[[r]] <- z[[r]] %*% Eb[1, 1:qri[r]]
        if (K_b > 1) for (k in 2:K_b) fit[[r]] <- cbind(fit[[r]], z[[r]] %*% Eb[k, 1:qri[r]])        
      }  

      if (pfi[r]){                           ### there are also fixed effects

        if (r == 1) alpha_r <- alpha[1:pfi[r]]
        else        alpha_r <- alpha[(sum(pfi[1:(r-1)])+1):sum(pfi[1:r])]
        
        if (object$fixed.intercept[r]){
          fit[[r]] <- fit[[r]] + alpha_r[1]
          alpha_r <- alpha_r[-1]
        }

        if (object$p[r]){
          if (object$p[r] == 1){
            if (is.matrix(x[[r]])) if (ncol(x[[r]]) != 1) stop("x[[", r, "]] must have 1 column", sep="")             
            x[[r]] <- matrix(x[[r]], ncol=1)
          }
          
          if (!is.matrix(x[[r]])) stop("x[[", r, "]] must be a matrix", sep="")
          if (ncol(x[[r]]) != object$p[r]) stop("x[[", r, "]] must have ", object$p[r], " columns", sep="")
          
          if (object$q[r]){      ## fit[[r]] is n x K matrix
            fit[[r]] <- fit[[r]] + matrix(rep(x[[r]] %*% alpha_r, K_b), ncol=K_b)
          }else{                 ## only random intercept among random effects -> fit[[r]] is a 1 x K matrix
            fit[[r]] <- matrix(rep(fit[[r]], nrow(x[[r]])), byrow=TRUE, ncol=K_b) + matrix(rep(x[[r]] %*% alpha_r, K_b), ncol=K_b)
          }  
        }  
      }

      ### inverse link function
      if (object$dist[r] == "binomial(logit)"){
        efit <- exp(fit[[r]])
        fit[[r]] <- efit / (1 + efit)
      }else{
        if (object$dist[r] == "poisson(log)"){
          fit[[r]] <- exp(fit[[r]])
        }else{
          if (object$dist[r] != "gaussian") stop("Not (yet) implemented for dist: ", object$dist[r], sep="")
        }  
      }  
      
      ##### glmer calculation if requested
      if (glmer){    ### if (glmer)

        if (object$dist[r] == "gaussian") next      ### Nothing must be done for Gaussian response.
                
        ### There is only random intercept in a model, nothing else otherwise.
        if (object$random.intercept[r] & object$q[r] == 0 & object$p[r] == 0){
          mfit[[r]] <- matrix(NA, nrow=1, ncol=K_b)                    
          for (k in 1:K_b){

            SDbk <- sqrt(as.numeric(Vb[[k]]))            
            
            if (object$dist[r] == "binomial(logit)"){

              ### Integrand
              Integrand <- function(b){
                return((exp(b) / (1 + exp(b))) * dnorm(b, mean = Eb[k,], sd = SDbk))
              }

              ### Mode of the integrand including the Hessian
              #ModeIntegr <- optimize(f = Integrand, interval = Eb[k,] + c(-7, 7)*SDbk, maximum = TRUE)
              #muHat <- ModeIntegr$maximum
              #sigmaHat <- SDbk
              
              ModeIntegr <- optim(par = Eb[k,], fn = Integrand, hessian = TRUE, control = list(fnscale = -1), method = "Brent", lower = Eb[k,] - 10*SDbk, upper = Eb[k,] + 10*SDbk)
              muHat <- ModeIntegr$par
              Hess <- as.numeric(ModeIntegr$hessian)
              if (Hess < 0) sigmaHat <- 1/sqrt(-Hess) else sigmaHat <- SDbk

              ### Integral by the Gauss-Hermite quadrature
              mfit[[r]][1, k] <- fastGHQuad::aghQuad(Integrand, muHat = muHat, sigmaHat = sigmaHat, rule = QRule)
              
            }else{
              if (object$dist[r] == "poisson(log)"){
                mfit[[r]][1, k] <- fit[[r]][1, k] * exp(as.numeric(Vb[[k]])/2)      ## = moment generating function of N in t=1
              }else{
                stop("glmer calculation [PART 2] not (yet) implemented for dist: ", object$dist[r], sep="")
              }  
            }  
          }            
        }else{

          ### Only random effects in a model (including random intercept)
          if (object$random.intercept[r] & object$q[r] > 0 & object$p[r] == 0){
            mfit[[r]] <- matrix(NA, nrow=nrow(z[[r]]), ncol=K_b)
            for (k in 1:K_b){                        
              if (object$dist[r] == "binomial(logit)"){
                
                for (i in 1:nrow(mfit[[r]])){

                  ### Integrand
                  Integrand <- function(b){
                    eta.random <- as.numeric(crossprod(z[[r]][i,], b))
                    return((exp(eta.random) / (1 + exp(eta.random))) * dMVN(b, mean = Eb[k,], Sigma = Vb[[k]]))
                  }

                  ### Mode of the integrand including the Hessian
                  ModeIntegr <- optim(par = Eb[k,], fn = Integrand, hessian = TRUE, control = list(fnscale = -1))
                  muHat <- ModeIntegr$par
                  Hess <- ModeIntegr$hessian

                  ### Integral by the Laplace approximation
                  stop("glmer calculation [random.intercept & q > 0 & p = 0] not yet implemented")                  
                  #mfit[[r]][i, k] <- 
                }
              }else{
                if (object$dist[r] == "poisson(log)"){
                  Z <- cbind(1, z[[r]])
                  LVb <- t(chol(Vb[[k]]))
                  ZLVb <- Z %*% LVb
                  ZVbZ <- apply(ZLVb, 1, crossprod)
                  mfit[[r]][, k] <- fit[[r]][, k] * exp(ZVbZ/2)      ## = moment generating function of N in t=(1, z)
                }else{
                  stop("glmer calculation [PART 3] not (yet) implemented for dist: ", object$dist[r], sep="")
                }  
              }  
            }
          }else{
            ### Random intercept is the only random effect in a model + there are some fixed effects as well
            if (object$random.intercept[r] & object$q[r] == 0 & object$p[r] > 0){
              mfit[[r]] <- matrix(NA, nrow=nrow(x[[r]]), ncol=K_b)
              
              for (k in 1:K_b){                        
                if (object$dist[r] == "binomial(logit)"){

                  SDbk <- sqrt(as.numeric(Vb[[k]]))
                  
                  for (i in 1:nrow(mfit[[r]])){

                    ### Integrand
                    eta.fixed <- as.numeric(crossprod(x[[r]][i,], alpha_r))
                    Integrand <- function(b){
                      return((exp(b + eta.fixed) / (1 + exp(b + eta.fixed))) * dnorm(b, mean = Eb[k,], sd = SDbk))
                    }

                    ### Mode of the integrand including the Hessian
                    #ModeIntegr <- optimize(f = Integrand, interval = Eb[k,] + c(-10, 10)*SDbk, maximum = TRUE)
                    #muHat <- ModeIntegr$maximum
                    #sigmaHat <- SDbk
                
                    ModeIntegr <- optim(par = Eb[k,], fn = Integrand, hessian = TRUE, control = list(fnscale = -1), method = "Brent", lower = Eb[k,] - 10*SDbk, upper = Eb[k,] + 10*SDbk)
                    muHat <- ModeIntegr$par
                    Hess <- as.numeric(ModeIntegr$hessian)
                    if (Hess < 0) sigmaHat <- 1/sqrt(-Hess) else sigmaHat <- SDbk

                    ### Integral by the Gauss-Hermite quadrature
                    mfit[[r]][i, k] <- fastGHQuad::aghQuad(Integrand, muHat = muHat, sigmaHat = sigmaHat, rule = QRule)
                  }
                  
                }else{
                  if (object$dist[r] == "poisson(log)"){
                    mfit[[r]][, k] <- fit[[r]][, k] * exp(as.numeric(Vb[[k]])/2)      ## = moment generating function of N in t=1
                  }else{
                    stop("glmer calculation [PART 3] not (yet) implemented for dist: ", object$dist[r], sep="")
                  }  
                }  
              }              
            }else{              
              if (object$random.intercept[r] & object$q[r] > 0 & object$p[r] > 0){
                stop("glmer calculation [random.intercept & q > 0 & p > 0] not yet implemented")
              }else{
                if (!object$random.intercept[r] & object$q[r] > 0 & object$p[r] == 0){
                  stop("glmer calculation [!random.intercept & q > 0 & p == 0] not yet implemented")
                }else{
                  if (!object$random.intercept[r] & object$q[r] > 0 & object$p[r] > 0){
                    stop("glmer calculation [!random.intercept & q > 0 & p > 0] not yet implemented")
                  }else{
                    stop("Programming error, please, contact AK.")
                  }  
                }  
              }  
            }              
          }  
        }

        fit[[r]] <- mfit[[r]]
      }              ### end of if (glmer)
      
    ### only fixed effects in the model  
    }else{
      stop("This part not yet implemented. Please, contact the author.")
    }

  }     ### end of loop (rr)
  
  return(fit)  
}  

