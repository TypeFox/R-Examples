## Functions to generate the models used in MRSP.
## Wolfgang Pößnecker, last modified: 10.11.2014, 23:19

MRSP.model <- function(invlink, link, logl, loglik, gradient, fisher,
                       constraint, check, sancheck, name = "user-specified",
                       comment = "user-specified"){
  ## Purpose: Generates models to be used for the MRSP package.
  ## ----------------------------------------------------------------------
  ## Arguments: in general, the functions take a generic data argument dat with
  ##            slots depending on the particular model and its data structure.
  ##            for MRSP, we need a list with elements "y", "x" and "V"
  ## 
  ## invlink:   a function with arguments "eta" implementing the inverse link
  ##            function. eta must be a matrix of size nobs * K for K response
  ##            categories.
  ## link:      a function with arguments "mu" and "constr" implementing the link
  ##            function for the various options of the constraint in "constr". 
  ## logl:      a function with arguments "y" and "mu" implementing
  ##            the log-likelihood contribution of each observation.  weights
  ##            must be a vector of length nobs. output is nobs x 1.
  ## loglik:    a function with arguments "y", "mu" and "weights"
  ##            implementing the log-likelihood function. weights must be a
  ##            vector of length nobs.
  ## gradient:  a function with arguments "x", "V", "y", "mu", "coef", "weights",
  ##            "penweights", "grpindex" and "Proximal.args" implementing the
  ##            gradient of the log-likelihood function. returns a list with two
  ##            elements, the first one being a K x p matrix that belongs to the
  ##            global predictors in x and the second a K x L matrix belonging
  ##            to the L category-specific predictors in V. x must be an nobs*p
  ##            matrix, V a list of length K with entries that are nobs x L
  ##            matrices. y must be a nobs x K matrix.
  ## fisher :   a function with arguments "x", "V", "mu" and "weights"
  ##            implementing the  hessian of the log-likelihood. same arguments
  ##            as above, returns a fisher matrix of size (K*p + K*L) x (K*p + L)
  ## constraint: the default constraint to use with the respective model class.
  ## sancheck:  a function with arguments "mu", "eta" and "weights" that checks
  ##            whether the fit converges towards a degenerated 0-1-separation. 
  ## name:      a character name

  RET <- new("MRSP.model",
           invlink    = invlink,
           link       = link,
           logl       = logl,
           loglik     = loglik,
           gradient   = gradient,
           fisher     = fisher,
           constraint = constraint,
           name       = name,
           check      = check,
           sancheck   = sancheck,
           comment    = comment)
   RET
}



## the model used for this package:
multinomlogit <- function(){
  MRSP.model(invlink  = cmpfun(function(eta){
                         m <- maxRow(eta)
                         eeta <- exp(eta - m)
                         mu <- eeta / rowSums(eeta)
                         return(mu)
                        }),
             link     = cmpfun(function(mu, constr){
                         if(constr == "none"){
                          eta <- log(mu)
                         }else{
                          if(is.numeric(constr)){
                           mu.ref <- mu[,constr, drop=T]
                          }else if(constr == "symmetric"){
                           mu.ref <- apply(mu, 1, function(u) exp(mean(log(u))))
                          }    
                          eta <- log(mu / mu.ref)
                         }
                         return(eta)
                        }),
             logl     = cmpfun(function(y, mu, ...){
                         mu[which(mu <= 1e-6)] <- 1e-8
                         mu[which(mu >= 1 - 1e-6)] <- 1 - 1e-8
                         out <- rowSums(y * log(mu))
                         return(out)
                        }),
             loglik   = cmpfun(function(y, mu, weights, ...){
                         mu[which(mu <= 1e-6)] <- 1e-8
                         mu[which(mu >= 1 - 1e-6)] <- 1 - 1e-8
                         out <- sum(weights * y * log(mu))
                         return(out)
                        }),
             gradient = cmpfun(function(dat, mu, coef, weights, penweights,
                                        grpindex, Proximal.args, ...){
                         constr <- Proximal.args$constraint
                         indg <- Proximal.args$indg
                         indcs <- Proximal.args$indcs
                         grad <- list()
                         grad[[1]] <- crossprod((dat$y - mu), weights * dat$x)
                         if(is.numeric(constr)) grad[[1]][constr,] <- 0
                         ## maybe mean-centering is necessary for constr 'symmetric', but probably not...
                         names(grad)[[1]] <- "global"
                         if(Proximal.args$hasV){

                          grad[[2]] <- matrix(nrow=ncol(dat$y), ncol=ncol(dat$V[[1]]))
                          if(length(indg) > 0){
                           gradg <- crossprod(c(dat$y - mu), rep(weights, times=ncol(dat$y)) * do.call("rbind", dat$V))[,indg]
                           grad[[2]][,indg] <- matrix(rep(gradg, ncol(dat$y)), nrow=ncol(dat$y), byrow=T)
                          }
                          if(length(indcs > 0)){
                           grad[[2]][,indcs] <- (t(mapply(crossprod, as.list(as.data.frame(dat$y - mu)), dat$V)))[,indcs]
                          }
                          names(grad)[[2]] <- "cat-specific"
                         }
                         return(grad)
                        }),
             fisher   = cmpfun(function(dat, mu, weights, Proximal.args, ...){
                         warning(paste("Fisher matrix needs rework and can't be used for the moment"))
                         #Fisher <- 0
                         #for(i in seq_len(nrow(dat$y))){
                         # sigmai <- diag(mu[i,]) - mu[i,] %*% t(mu[i,])
                         # Xi <- t(diag(ncol(mu)) %x% dat$x[i,])
                         # if(!is.null(dat$V)){Xi <- cbind(Xi, dat$V[[i]])}
                         # Fisher <- Fisher + crossprod(Xi, weights[i] * sigmai %*% Xi)
                         #}
                         #return(Fisher)
                        }),
             constraint = "symmetric",
             check    = function(dat){all(dat$y >= 0) & all(dat$y <= 1)},
             sancheck = function(coef, coef.old2, mu, eta, weights, Proximal.args){
                         if(Proximal.args$ridgestabil == T){Proximal.args$do.sancheck <<- F}
                         if(mean(abs(do.call("c", lapply(coef, "c")))) > mean(abs(do.call("c", lapply(coef.old2, "c"))))){
                          Proximal.args$sancount <<- Proximal.args$sancount + 1
                         }else{Proximal.args$sancount <<- 0}
                         if((mean(abs(do.call("c", lapply(coef, "c")))) > 5) & (Proximal.args$sancount > 5)){
                          if(!Proximal.args$ridgestabil){
                           warning(paste("The specified model seems to allow for a perfect separation of the response categories. A small ridge penalty was added to prevent diverging coefficients."))
                          } 
                          Proximal.args$ridgestabil <<- T
                          Proximal.args$ridgestabilrf <<- T
                          Proximal.args$do.sancheck <<- F
                         }
                        },
             name     = "Multinomial Logit Model",
             comment  = "Here y denotes the scaled response (i.e. y in [0,1]), if multiple observations of the same covariates are available, they are accounted for by the use of 'weights'"
  )
}





sequentiallogit <- function(){
             link     = cmpfun(function(mu, constr){
                         if(any(mu < -0.01)){
                          print(mu)
                          stop("mu negative in link")
                         } 
                         invmu <- 1 - rowCumsum(mu)
                         invmu <- cbind(1, invmu[,-ncol(invmu)])
                         invmu[which(abs(invmu) <= 1e-8)] <- 1e-8 
                         muinvmu <- mu / invmu     
                         muinvmu[which(muinvmu <= 1e-8)] <- 1e-8
                         muinvmu[which(muinvmu >= 1 - 1e-8)] <- 1 - 1e-8
                         eta <- log(muinvmu / (1 - muinvmu))
                         return(eta)
                        })
  MRSP.model(invlink  = cmpfun(function(eta){
                         lpeta <- log1pexp(eta)
                         if(any(is.na(lpeta))){
                          print(eta)
                          print(lpeta)
                          stop("lpeta is NA or NaN in invlink")
                         }
                         if(any(!is.finite(lpeta))){
                          print(lpeta)
                          stop("lpeta wrong in invlink")
                         } 
                         lmu <- eta - rowCumsum(lpeta)
                         if(any(lmu > 0.01)){
                          print(round(lmu, 4))
                          stop("'lmu' too large in invlink")
                         } 
                         mu <- exp(lmu)
                         return(mu)
                        }),
             link     =  link,
             logl     = cmpfun(function(y, mu, ...){
                         y <- cbind(y, 1 - rowSums(y))
                         mu <- cbind(mu, 1 - rowSums(mu))
                         mu[which(mu <= 1e-8)] <- 1e-8
                         mu[which(mu >= 1 - 1e-8)] <- 1 - 1e-8
                         out <- rowSums(y * log(mu))
                         return(out)
                        }),
             loglik   = cmpfun(function(y, mu, weights, ...){
                         y <- cbind(y, 1 - rowSums(y))
                         mu <- cbind(mu, 1 - rowSums(mu))
                         mu[which(mu <= 1e-8)] <- 1e-8
                         mu[which(mu >= 1 - 1e-8)] <- 1 - 1e-8
                         out <- sum(weights * y * log(mu))
                         return(out)
                        }),
             gradient = cmpfun(function(dat, mu, coef, weights, penweights,
                                        grpindex, Proximal.args, ...){
                         ## first the 'traditional' gradient/score function:
                         sl.indg <- Proximal.args$sl.indg
                         sl.indcs <- Proximal.args$sl.indcs
                         eta <- link(mu)
                         eta[which(eta >= 700)] <- 700
                         eeta <- exp(eta)
                         eeta1 <- 1 + eeta
                         maxy <- max.col(cbind(dat$y, 1 - rowSums(dat$y)))
                         grad <- list()
                         grad[[1]] <- matrix(nrow=ncol(dat$y), ncol=ncol(dat$x))
                         if(any(!is.finite(eeta))){
                          print(sl.indcs)
                          print(sl.indg)
                          stop("eeta not finite in gradient")
                         }
                         if(length(sl.indcs) > 0){
                          ytildecs <- yseqlog.constructor(-eeta, eeta1, maxy-1)
                          grad[[1]][,sl.indcs] <- crossprod(ytildecs, weights * dat$x[,sl.indcs])
                         } 
                         if(length(sl.indg) > 0){
                          etafrac <- eeta / eeta1
                          if(any(!is.finite(etafrac))){
                           print(etafrac)
                           stop("etafrac not finite in gradient for global coefs")
                          }
                          ceeta <- rowCumsum(etafrac)
                          ytildeg <- cbind(1 - ceeta, -ceeta[,ncol(dat$y)])
                          ytildeg <- ytildeg[cbind(seq_len(nrow(ytildeg)), maxy)]
                          gradg <- crossprod(ytildeg, weights * dat$x[,sl.indg])
                          grad[[1]][,sl.indg] <- matrix(rep(gradg, ncol(dat$y)), nrow=ncol(dat$y), byrow=T)
                         }
                         return(grad)
                        }),
             fisher   = cmpfun(function(dat, mu, weights, Proximal.args, ...){
                         warning(paste("Fisher matrix needs rework and can't be used for the moment"))
                        }),
             constraint = "none",
             check    = function(dat){all(dat$y >= 0) & all(dat$y <= 1)},
             sancheck = function(coef, coef.old2, mu, eta, weights, Proximal.args){
                         if(Proximal.args$ridgestabil == T){Proximal.args$do.sancheck <<- F}
                         if(mean(abs(do.call("c", lapply(coef, "c")))) > mean(abs(do.call("c", lapply(coef.old2, "c"))))){
                          Proximal.args$sancount <<- Proximal.args$sancount + 1
                         }else{Proximal.args$sancount <<- 0}
                         if((mean(abs(do.call("c", lapply(coef, "c")))) > 5) & (Proximal.args$sancount > 5)){
                          if(!Proximal.args$ridgestabil){
                           warning(paste("The specified model seems to allow for a perfect separation of the response categories. A small ridge penalty was added to prevent diverging coefficients."))
                          } 
                          Proximal.args$ridgestabil <<- T
                          Proximal.args$ridgestabilrf <<- T
                          Proximal.args$do.sancheck <<- F
                         }
                        },
             name     = "Sequential Logit Model",
             comment  = "Here y denotes the scaled response (i.e. y in [0,1]), if multiple observations of the same covariates are available, they are accounted for by the use of 'weights'")
}




CUBbinomiallogit <- function(){
  MRSP.model(invlink  = cmpfun(function(eta){
                         #eeta <- exp(-eta)
                         #mu <- 1 / (1 + eeta)
                         logeta <- eta - log1pexp(eta) # see MRSP-helpers for the definition of log1pexp
                         mu <- exp(logeta)
                         return(mu)
                        }),
             link     = cmpfun(function(mu, constr="none"){
                         if(constr == "none"){
                          eta <- log(mu/(1-mu))
                         }else{
                          stop("constraint must be 'none' for CUBbinomiallogit models")
                         }
                         return(eta)
                        }),
             logl     = cmpfun(function(y, mu, ...){
                         k <- max(y)
                         binomcoefs <- log(choose(k-1, y-1))
                         mu[which(mu <= 1e-6)] <- 1e-8
                         mu[which(mu >= 1 - 1e-6)] <- 1 - 1e-8
                         out <- rowSums( (y-1)*log(1-mu) + (k-y)*log(mu) + binomcoefs )
                         return(out)
                        }),
             loglik   = cmpfun(function(y, mu, weights, ...){
                         k <- max(y)
                         binomcoefs <- log(choose(k-1, y-1))
                         mu[which(mu <= 1e-6)] <- 1e-8
                         mu[which(mu >= 1 - 1e-6)] <- 1 - 1e-8
                         out <- sum(weights*( (y-1)*log(1-mu) + (k-y)*log(mu) + binomcoefs ))
                         return(out)
                        }),
             ## compute the score function. output must be a list for technical reasons.
             ## it is typically of length 1, with the entry being a K x p matrix
             ## in models with p-dimensional predictor dat$x and K-dimensional response
             ## variable dat$y.
             gradient = cmpfun(function(dat, mu, coef, weights, penweights,
                                        grpindex, Proximal.args, ...){
                         constr <- Proximal.args$constraint
                         k <- max(dat$y)
                         grad <- list()
                         grad[[1]] <- crossprod(((k-dat$y) - (k-1)*mu), weights * dat$x)
                         if(constr != "none") stop("constraint must be 'none' for CUBbinomiallogit models")
                         ## maybe mean-centering is necessary for constr 'symmetric', but probably not...
                         names(grad)[[1]] <- "global"
                         if(Proximal.args$hasV){
                          stop("function 'gradient' does currently not support 'V-type' predictors in family CUBbinomiallogit")
                         }
                         return(grad)
                        }),
             fisher   = cmpfun(function(dat, mu, weights, Proximal.args, ...){
                         warning(paste("Fisher matrix needs rework and can't be used at the moment"))
                        }),
             constraint = "none",
             ## the following function checks if all y-values are integers lying between 1 and k.
             check    = function(dat){k <- max(dat$y); all(dat$y %in% seq(from=1, to=k)) },#&& length(table(dat$y)) == k},
             ## NEVER EVER touch the following sanity check function which is used to detect and handle diverging coefficients.
             sancheck = function(coef, coef.old2, mu, eta, weights, Proximal.args){
                         if(Proximal.args$ridgestabil == T){Proximal.args$do.sancheck <<- F}
                         if(mean(abs(do.call("c", lapply(coef, "c")))) > mean(abs(do.call("c", lapply(coef.old2, "c"))))){
                          Proximal.args$sancount <<- Proximal.args$sancount + 1
                         }else{Proximal.args$sancount <<- 0}
                         if((mean(abs(do.call("c", lapply(coef, "c")))) > 5) & (Proximal.args$sancount > 5)){
                          if(!Proximal.args$ridgestabil){
                           warning(paste("The specified model seems to yield divering parameter estimates. A small ridge penalty was added to prevent diverging coefficients."))
                          }
                          Proximal.args$ridgestabil <<- T
                          Proximal.args$ridgestabilrf <<- T
                          Proximal.args$do.sancheck <<- F
                         }
                        },
             name     = "CUB Binomial Logit Model",
             comment  = "Here y takes integer values between 1 and k := max(y)"
  )
}




## the model used for this package:
OLSreg <- function(){
  MRSP.model(invlink  = cmpfun(function(eta){
                         mu <- eta
                         return(mu)
                        }),
             link     = cmpfun(function(mu, constr){
                         eta <- mu
                         return(eta)
                        }),
             logl     = cmpfun(function(y, mu, ...){
                         out <- -(y - mu)^2
                         return(out)
                        }),
             loglik   = cmpfun(function(y, mu, weights, ...){
                         out <- -sum(weights * (y - mu)^2)
                         return(out)
                        }),
             gradient = cmpfun(function(dat, mu, coef, weights, penweights,
                                        grpindex, Proximal.args, ...){
                         grad <- list()
                         grad[[1]] <- crossprod((dat$y - mu), weights * dat$x)
                         ## maybe mean-centering is necessary for constr 'symmetric', but probably not...
                         names(grad)[[1]] <- "global"
                         return(grad)
                        }),
             fisher   = cmpfun(function(dat, mu, weights, Proximal.args, ...){
                         warning(paste("Fisher matrix needs rework and can't be used for the moment"))
                         #Fisher <- 0
                         #for(i in seq_len(nrow(dat$y))){
                         # sigmai <- diag(mu[i,]) - mu[i,] %*% t(mu[i,])
                         # Xi <- t(diag(ncol(mu)) %x% dat$x[i,])
                         # if(!is.null(dat$V)){Xi <- cbind(Xi, dat$V[[i]])}
                         # Fisher <- Fisher + crossprod(Xi, weights[i] * sigmai %*% Xi)
                         #}
                         #return(Fisher)
                        }),
             constraint = "none",
             check    = function(dat){is.numeric(dat$y)},
             sancheck = function(coef, coef.old2, mu, eta, weights, Proximal.args){
                         if(Proximal.args$ridgestabil == T){Proximal.args$do.sancheck <<- F}
                         if(mean(abs(do.call("c", lapply(coef, "c")))) > mean(abs(do.call("c", lapply(coef.old2, "c"))))){
                          Proximal.args$sancount <<- Proximal.args$sancount + 1
                         }else{Proximal.args$sancount <<- 0}
                         if((mean(abs(do.call("c", lapply(coef, "c")))) > 5) & (Proximal.args$sancount > 5)){
                          if(!Proximal.args$ridgestabil){
                           warning(paste("The specified model seems to allow for a perfect separation of the response categories. A small ridge penalty was added to prevent diverging coefficients."))
                          }
                          Proximal.args$ridgestabil <<- T
                          Proximal.args$ridgestabilrf <<- T
                          Proximal.args$do.sancheck <<- F
                         }
                        },
             name     = "OLS Regression",
             comment  = "Linear OLS Regression. Note that functions 'logl' and 'loglik' are based on the RSS, not the actual loglikelihood! Also note that error variance estimation in penalized models can be subject to substantial bias."
  )
}