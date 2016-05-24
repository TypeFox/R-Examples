# Author: Harry Southworth
# Date: 2011-11-25
## Purpose: Create a predict method for objects of class evmOpt, evmSim
##          and evmBoot that
##          returns parameters, return levels or (maybe) return periods,
##          depending on arguments given.
#
# predict.evm
# predict.evmSim
# predict.evmBoot
# rl
# rl.evm
# rl.evmSim
# rl.evmBoot
# linearPredictors
# linearPredictors.evm
# linearPredictors.evmSim
# linearPredictors.evmBoot

################################################################################
## evm

predict.evmOpt <-
    # Get predictions for an evm object. These can either be the linear predictors
    # or return levels.
function(object, M=1000, newdata=NULL, type="return level", se.fit=FALSE,
         ci.fit=FALSE, alpha=.050, unique.=TRUE, ...){
    theCall <- match.call()

    res <- switch(type,
                  "rl"=, "return level" = rl.evmOpt(object, M, newdata,
                                                 se.fit=se.fit, ci.fit=ci.fit,
                                                 alpha=alpha, unique.=unique.),
                  "lp" =,"link" = linearPredictors.evmOpt(object, newdata, se.fit,
                                                   ci.fit, alpha, unique.=unique.)
                  )
    res
}

## Linear predictor functions for GPD

linearPredictors.evmOpt <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                             alpha=.050, unique.=TRUE, full.cov=FALSE, ...){

    D <- texmexMakeNewdataD(object, newdata)

    if (unique.){
        z <- do.call('cbind', D)
        u <- !duplicated(z)
        D <- lapply(D, function(x, u) {
                           if(is.matrix(x[u,]))  x[u, ]
                           else if(ncol(x) == 1)  cbind(x[u,])
                           else t(cbind(x[u,]))
                       }, u=u )
    }

    res <- texmexMakeParams(coef(object), D)
    colnames(res) <- names(D)

    # Get the covariance matrices - one for every unique observation
    if(ci.fit | se.fit | full.cov){
      cov.se <- texmexMakeCovariance(object$cov, D)
      # Get standard errors
      ses <- t(sapply(cov.se, function(x){ sqrt(diag(x)) }))
      colnames(ses) <- paste(colnames(res), '.se', sep = '')
    }

    if (ci.fit){
        ci <- texmexMakeCI(res, ses, alpha)
        res <- cbind(res, ci)
    } # Close if(ci.fit

    if (se.fit){
        res <- cbind(res, ses)
    } # Close if(se.fit

    for (i in 1:length(D)){
      res <- addCov(res, D[[i]])
    }

    res <- list(link=res,family=object$family)
    
    if (full.cov){
        res$cov <- cov.se
    }

    oldClass(res) <- "lp.evmOpt"
    res
}

## Return level functions for GPD

## Reversing arguments M and newdata for anyone who wants to call these functions
## directly

## Will want to get return levels when using GEV rather than GPD, so make
## rl generic

rl <- function(object, M = 1000, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, ...){
    UseMethod("rl")
}

linearPredictors <- function(object, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, ...){
    UseMethod("linearPredictors")
}


rl.evmOpt <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                       alpha=.050, unique.=TRUE, ...){
    co <- linearPredictors.evmOpt(object, newdata=newdata, unique.=unique., full.cov=TRUE)
    covs <- co$cov # list of covariance matrices, one for each (unique) observation
    co <- co$link
    X <- co[,-(1:length(object$data$D)), drop=FALSE]
    if(is.null(dim(X))){
      X <- matrix(X)
      dimnames(X) <- list(dimnames(co)[[1]], dimnames(co)[[2]][-(1:length(object$data$D))])
    }

    delta <- object$family$delta
    rl <- object$family$rl

    res <- lapply(M, rl, param=co, model=object)

    getse <- function(o, co, M, delta, covs){
        dxm <- lapply(split(co, 1:nrow(co)), delta, m=M, model=o)

        # Get (4.15) of Coles, page 82, adjusted for phi = log(sigma)
        se <- sapply(1:length(covs),
                     function(i, dxm, covs){
                        covs <- covs[[i]]; dxm <- c(dxm[[i]])
                        sqrt(mahalanobis(dxm, center=rep(0, ncol(covs)), cov=covs, inverted=TRUE))
                     }, dxm=dxm, covs=covs)
        se
    }

    if (ci.fit){
        ci.fun <- function(i, object, co, M, res, alpha, delta, covs){
            wh <- res[[i]];
            se <- getse(object, co, M[i], delta=delta, covs=covs)
            lo <- wh - qnorm(1 - alpha/2)*se
            hi <- wh + qnorm(1 - alpha/2)*se
            wh <- cbind(wh, lo=lo, hi=hi)

            colnames(wh) <- c("RL", paste(100*alpha/2, "%", sep = ""),
                              paste(100*(1 - alpha/2), "%", sep = ""))
            wh
        } # ci.fun
        res <- lapply(1:length(M), ci.fun, object=object, co=co,
                                           M=M, res=res, alpha=alpha,
                                           delta=delta, covs=covs)
    } # Close if (ci.fit

    if (se.fit){
        se.fun <- function(i, object, co, M, res, alpha, delta, covs){
            wh <- res[[i]]
            se <- getse(object, co, M[i], delta=delta, covs=covs)
            wh <- cbind(RL=wh, se.fit=se)
            wh
        } # ci.fun
        res <- lapply(1:length(M), se.fun, object=object, co=co,
                                           M=M, res=res, alpha=alpha,
                                           delta=delta, covs=covs)
    }

    cov.fun <- function(i,res){
      wh <- res[[i]]
      wh <- addCov(wh,X)
      wh
    }
    res <- lapply(1:length(M), cov.fun,res=res)

    names(res) <- paste("M.", M, sep = "")
    oldClass(res) <- "rl.evmOpt"
    res
}

################################################################################
## evmSim

predict.evmSim <- function(object, M=1000, newdata=NULL, type="return level",
                         se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE,
                         all=FALSE, sumfun=NULL, ...){
    theCall <- match.call()

    res <- switch(type,
                  "rl" = , "return level" = rl.evmSim(object, M=M, newdata=newdata,
                                                    se.fit=se.fit, ci.fit=ci.fit,
                                                    alpha=alpha, unique.=unique., all=all,
                                                    sumfun=sumfun,...),
                  "lp" = , "link" = linearPredictors.evmSim(object, newdata=newdata,
                                                      se.fit=se.fit, ci.fit=ci.fit,
                                                      alpha=alpha, unique.=unique., all=all,
                                                      sumfun=sumfun,...)
                  )
    res
}

linearPredictors.evmSim <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                                     alpha=.050, unique.=TRUE, all=FALSE, sumfun=NULL, ...){
    if (se.fit){ warning("se.fit not implemented - ignoring") }

    D <- texmexMakeNewdataD(object$map, newdata)

    X.all <- do.call("cbind", D)
    ModelHasCovs <- ncol(X.all) > length(D)

    if(ModelHasCovs){
      covCols <- apply(X.all, 2, function(x) !all(x==1))
      Xnames <- colnames(X.all)
      if(sum(covCols) == 1){
        X.all <- X.all[, covCols, drop=FALSE]
      }
      else {
        X.all <- X.all[, covCols]
      }
      colnames(X.all) <- Xnames[covCols]
    }

    if (unique.){
        X.all <- unique(X.all)
        u <- !duplicated(do.call("cbind", D))
        D <- lapply(D, function(x, u){ x[u,, drop=FALSE] }, u=u)
    }

    # Get matrices of parameters (i.e. split full parameter matrix into phi, xi whatever)
    param <- texmexGetParam(D, object$param)

    # Get linear predictors
    res <- lapply(1:nrow(D[[1]]), # For each observation get matrix of parameters
              function(i, x, p){
                  wh <- lapply(1:length(D),
                               function(j, x, p, i){
                                   rowSums(t(t(p[[j]]) * c(x[[j]][i, ])))
                               }, x=x, p=p, i=i)
                  wh <- do.call("cbind", wh)
                  colnames(wh) <- names(x)
                  wh
                }, x=D, p=param)
    # res should be a list containing a matrix for each observation.
    # The matrix represents the simulated posterior, one column for each
    # major parameter (i.e. linear predictors)

    ############################################################################
    ## Hard part should be done now. Just need to summarize

    if (ci.fit){
        # Need to get names by pasting together CI names and parameter names
        wh <- texmexMakeCISim(res[[1]], alpha=alpha, object=object, sumfun=sumfun)
        wh <- colnames(wh)
        wh <- paste(rep(names(D), ea=length(wh)), wh, sep = ": ")

        # Need to get order of output correct, so need to faff about with transposing
        res <- sapply(res, function(x){
                               t(texmexMakeCISim(x, alpha=alpha, object=object, sumfun=sumfun))
                           })
        res <- t(res)
        colnames(res) <- wh
    }

    else if (all){ res <- res }
    else { # Just point estimates
        res <- t(sapply(res, function(x){ apply(x, 2, mean) }))
    }
    if(!all){
      if(ModelHasCovs){
        for (i in 1:length(D)){
            res <- addCov(res,D[[i]])
        }
      }
    }
    else {
        if (ModelHasCovs & nrow(X.all) != length(res)){
            stop("Number of unique combinations of covariates doesn't match the number of parameters")
        }
        for (i in 1:length(res)){
          if(ModelHasCovs){
            res[[i]] <- cbind(res[[i]], matrix(rep(X.all[i,], nrow(res[[i]])),
                                               nrow=nrow(res[[i]]), byrow=TRUE))
            colnames(res[[i]]) <- c(names(D), colnames(X.all))
          } else {
            colnames(res[[i]]) <- names(D)
          }
        }
    }

    res <- list(link=res,family=object$map$family)
    oldClass(res) <- "lp.evmSim"
    res
}

rl.evmSim <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE, all=FALSE, sumfun=NULL,...){
    if (se.fit){ warning("se.fit not implemented") }

    co <- linearPredictors.evmSim(object, newdata=newdata, unique.=unique., all=TRUE, sumfun=NULL)$link
    # XXX Next line seems silly! Why not compute it from the line above?
    Covs <- linearPredictors.evmSim(object, newdata=newdata, unique.=unique., sumfun=NULL)$link
    X <- Covs[,-(1:length(object$map$data$D))]
    if(is.null(dim(X))){
      X <- matrix(X)
      dimnames(X) <- list(dimnames(Covs)[[1]],dimnames(Covs)[[2]][-(1:length(object$map$data$D))])
    }

    sim.rl <- function(m, param, model){
        rl <- model$family$rl
        rl(m=m, param, model)
    }

    # co is a list with one element for each unique item in
    # new data. Need to loop over vector M and the elements of co

    getrl <- function(m, co, ci.fit, alpha, all, object){
        res <- sapply(co, sim.rl, m=m, model=object$map)
        if (ci.fit){
            res <- texmexMakeCISim(res, alpha, object$map, sumfun, M=m)
        } # Close if (ci.fit
        else if (!all){
            res <- apply(res, 2, mean)
        }
        res
    }

    res <- lapply(M, getrl, co=co, ci.fit=ci.fit, alpha=alpha, all=all, object=object)

    if(!all){
      cov.fun <- function(i,res){
        wh <- res[[i]]
        wh <- addCov(wh,X)
        wh
      }
      res <- lapply(1:length(M), cov.fun,res=res)
    }

    names(res) <- paste("M.", M, sep = "")
    oldClass(res) <- "rl.evmSim"
    res
}

################################################################################
## evmBoot

predict.evmBoot <- function(object, M=1000, newdata=NULL, type="return level",
                            se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE,
                            all=FALSE, sumfun=NULL, ...){
    theCall <- match.call()

    res <- switch(type,
                  "rl" = , "return level" = rl.evmBoot(object, newdata=newdata, M=M,
                                                       se.fit=se.fit, ci.fit=ci.fit,
                                                       alpha=alpha, unique.=TRUE,
                                                       all=all, sumfun=sumfun,...),
                  "lp" = , "link" = linearPredictors.evmBoot(object, newdata=newdata,
                                                         se.fit=se.fit, ci.fit=ci.fit,
                                                         alpha=alpha, unique.=TRUE,
                                                         all=all, sumfun=sumfun,...)
                  )
    res
}

namesBoot2sim <- function(bootobject){
    names(bootobject) <- c("call", "param", "map")
    bootobject
}

linearPredictors.evmBoot <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=.050,
                                 unique.=TRUE, all=FALSE, sumfun=NULL,...){
    # This should just be the same as for an evmSim object, but some
    # names and stuff are different.
  object <- namesBoot2sim(object)
  res <- linearPredictors.evmSim(object, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit, all=all, unique.=unique., alpha=alpha, sumfun=sumfun,...)
  oldClass(res) <- "lp.evmBoot"
  res
}

rl.evmBoot <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=0.050, unique.=TRUE, all=FALSE, sumfun=NULL,...){
    # This should just be the same as for an evmSim object, but some
    # names are different.
  object <- namesBoot2sim(object)
  res <- rl.evmSim(object, M=M, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit,alpha=alpha, unique.=unique., all=all, sumfun=sumfun,...)
  oldClass(res) <- "rl.evmBoot"
  res
}

################################################################################
## Method functions

print.rl.evmOpt <- function(x, digits=3, ...){
    nms <- names(x)
    newnms <- paste("M =", substring(nms, 3), "predicted return level:\n")
    lapply(1:length(x), function(i, o, title){
                                 cat(title[i])
                                 print(o[[i]], digits=digits,...)
                                 cat("\n")
                                 NULL}, o=x, title=newnms)
    invisible(x)
}

summary.rl.evmOpt <- function(object, digits=3, ...){
    print.rl.evmOpt(object, digits=digits, ...)
}

print.rl.evmSim    <- print.rl.evmOpt
print.rl.evmBoot <- print.rl.evmOpt

summary.rl.evmSim    <- summary.rl.evmOpt
summary.rl.evmBoot <- summary.rl.evmOpt


print.lp.evmOpt <- function(x, digits=3, ...){
    cat("Linear predictors:\n")
    print(unclass(x$link), digits=3,...)
    invisible(x)
}

summary.lp.evmOpt <- function(object, digits=3, ...){
    print.lp.evmOpt(object, digits=3, ...)
}

summary.lp.evmSim    <- summary.lp.evmOpt
summary.lp.evmBoot <- summary.lp.evmOpt

print.lp.evmSim    <- print.lp.evmOpt
print.lp.evmBoot <- print.lp.evmOpt

