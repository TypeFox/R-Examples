# estimate
# inside control: lambda, lambdas, alpha.min, alpha.max, beta.max, iter, sd.thresh, scalex
bsgw <- function(formula, data, formulas=formula, weights, subset, na.action=na.fail, init="survreg"
                 , ordweib=FALSE, scale=0, control=bsgw.control(), print.level=2) {
  # TODO: implement weights, subset, na.action, scale
  # TODO: add robustness to init: 1) add bsgw objects, 2) allow failover when survreg doesn't work
  mycall <- match.call()
  if (!missing(weights)) warning("weights argument not supported yet, this argument will be ignored")
  if (!missing(subset)) warning("subset argument not supported yet, this argument will be ignored")
  if (!identical(na.action,na.fail)) stop("na.action argument not supported yet; only na.fail is currently accepted")
  if (scale>0) stop("scale argument not supported yet")
  
  if (ordweib) {
    formulas <- bsgw.strip.formula(formulas)
  }
  
  # TODO: need to make sure number of rows in X and Xs ends up being equal (watch out for removal of missing rows)
  mf <- model.frame(formula, data, drop.unused.levels=TRUE, na.action = na.fail) # incorporate na.action argument
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)
  y <- model.response(mf, "numeric")
  colnamesX <- colnames(X)
  if (colnamesX[1]!="(Intercept)") stop("intercept term must be included in formula")
  
  mfs <- model.frame(formulas, data, drop.unused.levels=TRUE, na.action = na.fail) # incorporate na.action argument
  mts <- attr(mfs, "terms")
  Xs <- model.matrix(mts, mfs)
  colnamesXs <- colnames(Xs)
  if (colnamesXs[1]!="(Intercept)") stop("intercept term must be included in formulas")

  if (init[1]=="survreg") { # using survreg to initialize the coefficients
    wreg <- survreg(formula, data) # TODO: make sure all relevant parameters are passed to this call, e.g. na.action
    alpha.ow <- 1/wreg$scale
    betas.0 <- if (ordweib) log(alpha.ow) else log((alpha.ow-control$alpha.min)/(control$alpha.max-alpha.ow))
    beta <- unname(-wreg$coefficients/wreg$scale)
    init <- list(beta=beta, betas=c(betas.0, rep(0, ncol(Xs)-1)))
    survreg.scale.ref <- wreg$scale
  } else if (class(init)=="survreg") { # using a previous survreg estimation object
    wreg <- init
    alpha.ow <- 1/wreg$scale
    betas.0 <- if (ordweib) log(alpha.ow) else log((alpha.ow-control$alpha.min)/(control$alpha.max-alpha.ow))
    beta <- unname(-wreg$coefficients/wreg$scale)
    init <- list(beta=beta, betas=c(betas.0, rep(0, ncol(Xs)-1)))
    survreg.scale.ref <- wreg$scale
  } else {
    warning("init argument not recognized, initializing all coefficient to zero") # TODO: is this desired behavior?
    init <- list(beta=rep(0,ncol(X)), betas=rep(0,ncol(Xs)))
    survreg.scale.ref <- NULL
    wreg <- NULL
  }
  
  if (is.list(control$scalex)) {
    X <- bsgw.scale(X, apply.sc=control$scalex$apply.scale.X, center=control$scalex$centerVec.X, scale=control$scalex$scaleVec.X)
    Xs <- bsgw.scale(Xs, apply.sc=control$scalex$apply.scale.Xs, center=control$scalex$centerVec.Xs, scale=control$scalex$scaleVec.Xs)
    control$scalex <- TRUE
  } else if (control$scalex) {
    X <- bsgw.scale(X)
    Xs <- bsgw.scale(Xs)
  }
  
  if (ordweib) formulas <- bsgw.strip.formula(formulas)
  
  ret <- list(call=mycall, formula=formula, formulas=formulas, weights=rep(1,nrow(X)), subset=1:nrow(X)
              , na.action=na.action, init=init, ordweib=ordweib, survreg.scale.ref=survreg.scale.ref, ordreg=wreg, scale=scale, control=control
              , X=X, Xs=Xs, y=y
              , contrasts=attr(X, "contrasts"), contrastss=attr(Xs, "contrasts")
              , xlevels=.getXlevels(mt, mf), xlevelss=.getXlevels(mts, mfs)
              , terms=mt, termss=mts
              , colnamesX=colnamesX, colnamesXs=colnamesXs
  )
  if (is.list(control$scalex) || control$scalex) {
    ret <- c(ret, list(apply.scale.X=attr(X, "apply.scale"), apply.scale.Xs=attr(Xs, "apply.scale")
                       , centerVec.X=attr(X, "centerVec"), scaleVec.X=attr(X, "scaleVec")
                       , centerVec.Xs=attr(Xs, "centerVec"), scaleVec.Xs=attr(Xs, "scaleVec")))
  }
  
  mcmc <- bsgw.mcmc(X, Xs, y[,1], y[,2], control$lambda, control$lambdas, iter=control$iter, sd.thresh=control$sd.thresh
                    , init=init, ordweib=ordweib, alpha.min=control$alpha.min, alpha.max=control$alpha.max, beta.max=control$beta.max
                    , print.level=print.level, nskip=control$nskip)
  sel <- (control$burnin+1):control$iter
  median <- list(beta=apply(mcmc$beta[sel,,drop=F], 2, median), betas=apply(mcmc$betas[sel,,drop=F], 2, median)
                 , survreg.scale=apply(mcmc$scale[sel,,drop=F], 2, median))
  
  km.fit <- survfit(bsgw.strip.formula(formula), data)
  
  ret <- c(ret, list(idx=mcmc$idx, idxs=mcmc$idx, median=median
              , smp=list(beta=mcmc$beta, betas=mcmc$betas, survreg.scale=mcmc$scale, lp=mcmc$lp, loglike=mcmc$loglike, logpost=mcmc$logpost)
           , km.fit=km.fit, tmax=max(y[,1])))
  class(ret) <- "bsgw"
  return (ret)
}

print.bsgw <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Scale formula:\n")
  print(x$formula)
  cat("Shape formula:\n")
  print(x$formulas)
  cat("lambda (shrinkage for scale coefficients):", x$control$lambda, "\n")
  cat("lambdas (shrinkage for shape coefficients):", x$control$lambdas, "\n")
  cat("ordinary (constant-shape) Weibull:", x$ordweib, "\n")
  if (!x$ordweib) {
    cat("Lower bound on shape parameter:", x$control$alpha.min, "\n")
    cat("Upper bound on shape parameter:", x$control$alpha.max, "\n")
  }
  cat("Number of MCMC iterations:", x$control$iter, "\n")
  cat("Number of burn-in iterations (used in calculating medians):", x$control$burnin, "\n")
  cat("Threshold on standard deviation of covariates:", x$control$sd.thresh, "\n")
  cat("Model matrix is scaled:", x$control$scalex, "\n")
  cat("Scale Coefficients:\n")
  beta <- x$median$beta; names(beta) <- x$colnamesX
  print(beta)
  cat("Shape Coefficients:\n")
  betas <- x$median$betas; names(betas) <- x$colnamesXs
  print(betas)
  cat("mean survreg-style scale parameter:", mean(x$median$survreg.scale), "\n")
  if (!is.null(x$survreg.scale.ref)) cat("\tscale parameter from non-Bayesian Weibull regression:", x$survreg.scale.ref, "\n")
  cat("number of observations:", nrow(x$X), "\n")
}

plot.bsgw <- function(x, pval=0.05, burnin=round(x$control$iter/2), nrow=2, ncol=3, ...) {
  iter <- x$control$iter
  sel <- (burnin+1):x$control$iter
  nsel <- length(sel)
  CI_prob <- c(pval/2, 0.5, 1-pval/2)
  nplot_per_page <- nrow*ncol
  # determine number of beta coefficients
  nbeta <- ncol(x$smp$beta)
  nbetas <- ncol(x$smp$betas)
  npage_beta <- ceiling(nbeta/nplot_per_page)
  npage_betas <- ceiling(nbetas/nplot_per_page)
  
  ## loglike and logpost
  # TODO: adjust ylim to make sure ordreg loglike is visible
  par(mfrow=c(2,2))
  plot(x$smp$loglike, type="l", xlab="Iteration", ylab="Log-likelihood", main="Log-likelihood, All")
  if (!is.null(x$ordreg)) abline(h=x$ordreg$loglik[2], lty=2, col="red")
  plot(x$smp$logpost, type="l", xlab="Iteration", ylab="Log-posterior", main="Log-posterior, All")
  if (!is.null(x$ordreg)) abline(h=x$ordreg$loglik[2], lty=2, col="red")
  plot(x$smp$loglike[sel], type="l", xlab="Iteration", ylab="Log-likelihood", main="Log-likelihood, Post-Burnin")
  if (!is.null(x$ordreg)) abline(h=x$ordreg$loglik[2], lty=2, col="red")
  plot(x$smp$logpost[sel], type="l", xlab="Iteration", ylab="Log-posterior", main="Log-posterior, Post-Burnin")
  if (!is.null(x$ordreg)) abline(h=x$ordreg$loglik[2], lty=2, col="red")
  
  # histogram of scale parameters for training set
  if (!x$ordweib) {
    par(mfrow=c(1,1))
    hist(x$median$survreg.scale, main="Survreg-Style Scale Parameter - Training Set", xlab="Survreg Scale Parameter")
  }
  
  ## traceplots
  # beta
  beta_q <- apply(x$smp$beta[sel,,drop=F], 2, quantile, probs=CI_prob)
  beta_lower <- beta_q[1,]
  beta_median <- beta_q[2,]
  beta_upper <- beta_q[3,]
  for (n in 1:npage_beta) {
    par(mfrow=c(nrow,ncol))
    offset <- (n-1)*nplot_per_page
    for (i in 1:nplot_per_page) {
      if (offset+i<=nbeta) {
        beta_ylim <- range(x$smp$beta[,offset+i], 0.0)
        plot(x$smp$beta[,offset+i], type="l", xlab="Iteration", ylab="Sample Value", ylim=beta_ylim
             , main = paste("beta[", x$colnamesX[offset+i], "]", sep=""))
        abline(h = 0)
        lines(sel, rep(beta_lower[offset+i], nsel), lty=2, col="red")
        lines(sel, rep(beta_median[offset+i], nsel), lty=2, col="red")
        lines(sel, rep(beta_upper[offset+i], nsel), lty=2, col="red")
      }
    }
  }

  if (x$ordweib) {
    if (!is.null(x$survreg.scale.ref)) yref <- x$survreg.scale.ref
    else yref <- 0.0
    survreg.scale_q <- quantile(x$smp$survreg.scale[sel,1], probs=CI_prob)
    survreg.scale_lower <- survreg.scale_q[1]
    survreg.scale_median <- survreg.scale_q[2]
    survreg.scale_upper <- survreg.scale_q[3]
    par(mfrow=c(1,1))
    survreg.scale_ylim <- range(x$smp$survreg.scale[,1], yref)
    plot(x$smp$survreg.scale[,1], type="l", xlab="Iteration", ylab="Sample Value", ylim=survreg.scale_ylim
         , main = paste("Survreg Scale", sep=""))
    abline(h = yref)
    lines(sel, rep(survreg.scale_lower, nsel), lty=2, col="red")
    lines(sel, rep(survreg.scale_median, nsel), lty=2, col="red")
    lines(sel, rep(survreg.scale_upper, nsel), lty=2, col="red")
  } else {
    # betas
    betas_q <- apply(x$smp$betas[sel,], 2, quantile, probs=CI_prob)
    betas_lower <- betas_q[1,]
    betas_median <- betas_q[2,]
    betas_upper <- betas_q[3,]
    for (n in 1:npage_betas) {
      par(mfrow=c(nrow,ncol))
      offset <- (n-1)*nplot_per_page
      for (i in 1:nplot_per_page) {
        if (offset+i<=nbetas) {
          betas_ylim <- range(x$smp$betas[,offset+i], 0.0)
          plot(x$smp$betas[,offset+i], type="l", xlab="Iteration", ylab="Sample Value", ylim=betas_ylim
               , main = paste("betas[", x$colnamesXs[offset+i], "]", sep=""))
          abline(h = 0)
          lines(sel, rep(betas_lower[offset+i], nsel), lty=2, col="red")
          lines(sel, rep(betas_median[offset+i], nsel), lty=2, col="red")
          lines(sel, rep(betas_upper[offset+i], nsel), lty=2, col="red")
        }
      }
    }
  }
  
  ## autocorrelation plots
  # beta
  for (n in 1:npage_beta) {
    par(mfrow=c(nrow,ncol))
    offset <- (n-1)*nplot_per_page
    for (i in 1:nplot_per_page) {
      if (offset+i<=nbeta) {
        if ((offset+i) %in% x$idx) {
          acf(x$smp$beta[sel,offset+i], main=paste("beta[", x$colnamesX[offset+i], "]", sep=""))
        } else {
          bsgw.empty.plot(main=paste("beta[", x$colnamesX[offset+i], "]", sep=""))
        }
      }
    }
  }
  
  if (x$ordweib) {
    par(mfrow=c(1,1))
    acf(x$smp$survreg.scale[sel], main=paste("Survreg Scale", sep=""))
  } else {
    # betas
    for (n in 1:npage_betas) {
      par(mfrow=c(nrow,ncol))
      offset <- (n-1)*nplot_per_page
      for (i in 1:nplot_per_page) {
        if (offset+i<=nbetas) {
          if ((offset+i) %in% x$idxs) {
            acf(x$smp$betas[sel,offset+i], main=paste("betas[", x$colnamesXs[offset+i], "]", sep=""))
          } else {
            bsgw.empty.plot(main=paste("betas[", x$colnamesXs[offset+i], "]", sep=""))
          }
        }
      }
    }
  }
  
  ## histograms
  # beta
  for (n in 1:npage_beta) {
    par(mfrow=c(nrow,ncol))
    offset <- (n-1)*nplot_per_page
    for (i in 1:nplot_per_page) {
      if (offset+i<=nbeta) {
        if ((offset+i) %in% x$idx) {
          hist(x$smp$beta[sel,offset+i], xlab="Sample Value"
               , main=paste("beta[", x$colnamesX[offset+i], "]", sep=""))
          abline(v = 0)
          abline(v = beta_median[offset+i], lty=2)
          abline(v = beta_lower[offset+i], lty=3)
          abline(v = beta_upper[offset+i], lty=3)
        } else {
          bsgw.empty.plot(main=paste("beta[", x$colnamesX[offset+i], "]", sep=""))
        }
      }
    }
  }
  
  # TODO: create branch for ordweib
  if (x$ordweib) {
    par(mfrow=c(1,1))
    hist(x$smp$survreg.scale[sel,1], xlab="Sample Value"
         , main=paste("Survreg Scale", sep=""))
    abline(v = 0)
    abline(v = survreg.scale_median, lty=2)
    abline(v = survreg.scale_lower, lty=3)
    abline(v = survreg.scale_upper, lty=3)
  } else {
    # betas
    for (n in 1:npage_betas) {
      par(mfrow=c(nrow,ncol))
      offset <- (n-1)*nplot_per_page
      for (i in 1:nplot_per_page) {
        if (offset+i<=nbetas) {
          if ((offset+i) %in% x$idxs) {
            hist(x$smp$betas[sel,offset+i], xlab="Sample Value"
                 , main=paste("betas[", x$colnamesXs[offset+i], "]", sep=""))
            abline(v = 0)
            abline(v = betas_median[offset+i], lty=2)
            abline(v = betas_lower[offset+i], lty=3)
            abline(v = betas_upper[offset+i], lty=3)
          } else {
            bsgw.empty.plot(main=paste("betas[", x$colnamesXs[offset+i], "]", sep=""))
          }
        }
      }
    }
  }
}

summary.bsgw <- function(object, pval=0.05, burnin=object$control$burnin, ...) {
  iter <- object$control$iter
  CI_prob <- c(pval/2, 0.5, 1-pval/2)
  sel <- (burnin+1):iter
  
  # beta
  beta_q <- apply(object$smp$beta[sel,,drop=F], 2, quantile, probs=CI_prob)
  beta_lb <- beta_q[1,]
  beta_med <- beta_q[2,]
  beta_ub <- beta_q[3,]
  beta_pval <- apply(object$smp$beta[sel,,drop=F], 2, bsgw.calc.pval, ref=0.0)
  coefficients_beta <- as.matrix(cbind(beta_med, beta_lb, beta_ub, beta_pval))
  dimnames(coefficients_beta) <- list(object$colnamesX, c("Estimate", "Lower Bound", "Upper Bound", "P-val"))
  # betas
  betas_q <- apply(object$smp$betas[sel,,drop=F], 2, quantile, probs=CI_prob)
  betas_lb <- betas_q[1,]
  betas_med <- betas_q[2,]
  betas_ub <- betas_q[3,]
  betas_pval <- apply(object$smp$betas[sel,,drop=F], 2, bsgw.calc.pval, ref=0.0)
  coefficients_betas <- as.matrix(cbind(betas_med, betas_lb, betas_ub, betas_pval))
  dimnames(coefficients_betas) <- list(object$colnamesXs, c("Estimate", "Lower Bound", "Upper Bound", "P-val"))
  
  # alpha (for training set)
  survreg.scale_q <- apply(object$smp$survreg.scale[sel,], 2, quantile, probs=CI_prob)
  survreg.scale_lb <- survreg.scale_q[1,]
  survreg.scale_med <- survreg.scale_q[2,]
  survreg.scale_ub <- survreg.scale_q[3,]
  
  ret <- list(call=object$call, pval=pval, burnin=burnin
              , coefficients=list(beta=coefficients_beta, betas=coefficients_betas)
              , survreg.scale=list(lower=survreg.scale_lb, median=survreg.scale_med, upper=survreg.scale_ub)
  )
  class(ret) <- "summary.bsgw"
  return (ret)
}

print.summary.bsgw <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("number of burn-in iterations discarded:", x$burnin, "\n")
  cat("confidence interval:", x$pval, "\n")
  cat("mean of median survreg-scale:", mean(x$survreg.scale$median), "\n")
  cat("## coefficients ##\n")
  cat("scale parameter:\n")
  print(x$coefficients$beta)
  cat("\nshape parameter:\n")
  print(x$coefficients$betas)
}

predict.bsgw <- function(object, newdata=NULL, tvec=NULL, burnin=object$control$burnin, ncores=1, ...) {
  iter <- object$control$iter
  alpha.min <- object$control$alpha.min
  alpha.max <- object$control$alpha.max
    
  tt <- object$terms
  tts <- object$termss
  Terms <- delete.response(tt)
  Termss <- delete.response(tts)
  
  if (is.null(newdata)) {
    nobs <- nrow(object$X)
    X <- object$X
    Xs <- object$Xs
    km.fit <- object$km.fit
  } else {
    newdata <- droplevels(newdata)
    mf <- model.frame(Terms, newdata, xlev = object$xlevels)
    mfs <- model.frame(Termss, newdata, xlev = object$xlevelss)
    X <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)
    Xs <- model.matrix(Termss, mfs, contrasts.arg = object$contrastss)
    
    if (object$control$scalex) {
      X <- bsgw.scale(X, apply.sc=object$apply.scale.X, center=object$centerVec.X, scale=object$scaleVec.X)
      Xs <- bsgw.scale(Xs, apply.sc=object$apply.scale.Xs, center=object$centerVec.Xs, scale=object$scaleVec.Xs)
    }
    nobs <- nrow(newdata)
  }

  if (!is.null(tvec)) {
  
    # TODO: we need an upper bound on length of tvec to avoid memory blow-up
    if (length(tvec)==1) tvec <- seq(from=0.0, to=object$tmax, length.out=tvec) # tvec is interpreted as number of time points
    nt <- length(tvec)
    mem.gb.per.iter <- 3*8*nt*nobs/(1024*1024*1024)
    required.mem.gb <- 3*8*iter*nt*nobs/(1024*1024*1024)
    
    tvec <- as.matrix(tvec)
    t_mat <- tvec[,rep(1,nobs)]
    
    # making sure returned prediction objects are not too big
    if (required.mem.gb>object$control$memlim.gb)
      stop("require memory exceeds specified limit\nconsider increasing limit via memlim.gb parameter or making prediction for a subset of data")
    
    if (ncores==1) {
      xbeta.all <- X%*%t(object$smp$beta)
      xbetas.all <- Xs%*%t(object$smp$betas)
      ret <- lapply(1:iter, FUN=function(i) {
        xbeta <- xbeta.all[,i,drop=F]
        xbetas <- xbetas.all[,i,drop=F]
        if (object$ordweib) {
          alpha <- exp(xbetas)
        } else {
          alpha <- alpha.min + (alpha.max-alpha.min)/(1+exp(-xbetas))
        }
        survreg.scale <- 1/alpha
        exbeta <- as.matrix(exp(xbeta))
        alpha_mat <- t(alpha[,rep(1,nt)])
        exbeta_mat <- t(exbeta[,rep(1,nt)])
        
        Htmp <- (t_mat^alpha_mat)*exbeta_mat
        htmp <- alpha_mat*(t_mat^(alpha_mat-1))*exbeta_mat
        
        return (list(h=htmp, H=Htmp, survreg.scale=survreg.scale))
      })
    } else {
      registerDoParallel(ncores)
      ret <- foreach (i=1:iter, .options.multicore=list(preschedule=TRUE)) %dopar% {
        xbeta <- X%*%object$smp$beta[i,]
        xbetas <- Xs%*%object$smp$betas[i,]
        if (object$ordweib) {
          alpha <- exp(xbetas)
        } else {
          alpha <- alpha.min + (alpha.max-alpha.min)/(1+exp(-xbetas))
        }
        survreg.scale <- 1/alpha
        exbeta <- as.matrix(exp(xbeta))
        alpha_mat <- t(alpha[,rep(1,nt)])
        exbeta_mat <- t(exbeta[,rep(1,nt)])
        
        Htmp <- (t_mat^alpha_mat)*exbeta_mat
        htmp <- alpha_mat*(t_mat^(alpha_mat-1))*exbeta_mat
        
        return (list(h=htmp, H=Htmp, survreg.scale=survreg.scale))
      }
    }

    h <- array(NA, dim=c(iter, nt, nobs))
    H <- array(NA, dim=c(iter, nt, nobs))
    survreg.scale <- array(NA, dim=c(iter, nobs))
    for (i in 1:iter) {
      h[i,,] <- ret[[i]]$h
      H[i,,] <- ret[[i]]$H
      survreg.scale[i,] <- ret[[i]]$survreg.scale
    }

    S <- exp(-H)
    
    survreg.scale_median <- apply(survreg.scale[(burnin+1):iter,], 2, median)
  
  } else {
    h <- NA
    H <- NA
    S <- NA
    survreg.scale <- NA
    survreg.scale_median <- NA
  }

  if (is.null(newdata)) {
    y <- object$y
    do_loglike <- T
  } else {
    Rterms <- drop.terms(tt)
    if (all(all.vars(Rterms)[1:2] %in% colnames(newdata))) {
      mfy <- model.frame(Rterms, newdata, xlev = object$xlevels) # TODO: add check to make sure response variable is available for newdata
      y <- model.response(mfy, "numeric")
      do_loglike <- T
      km.fit <- survfit(bsgw.strip.formula(object$formula), newdata)
    } else {
      do_loglike <- F
      km.fit <- NULL
    }
  }
  
  if (do_loglike) { # TODO: include logpost
    loglike <- sapply(1:iter, FUN=function(i) {
      bsgw.loglike(c(object$smp$beta[i,], object$smp$betas[i,]), X, Xs, y[,1], y[,2], object$ordweib, alpha.min, alpha.max)
    })
    loglike_median <- median(loglike[(burnin+1):iter])
  } else {
    loglike <- NA
    loglike_median <- NA
  }
  ret <- list(tvec=as.vector(tvec), burnin=burnin, median=list(loglike=loglike_median, survreg.scale=survreg.scale_median)
              , smp=list(h=h, H=H, S=S, survreg.scale=survreg.scale, loglike=loglike), km.fit=km.fit)
  class(ret) <- "predict.bsgw"
  return (ret)
}

summary.predict.bsgw <- function(object, idx=1:length(object$median$survreg.scale), burnin=object$burnin
                                 , pval=0.05, popmean=identical(idx,1:length(object$median$survreg.scale))
                                 , make.plot=TRUE, ...) { # TODO: verify this function, it was accidentally modified during editing of mixture model
  if (!all(idx %in% 1:length(object$median$survreg.scale))) {
    stop("invalid idx argument")
  }
  if (is.null(object$tvec)) {
    cat("prediction summary must be applied to time-dependent prediction entities")
    return (NULL) # TODO: can't we still return something useful?!
  }
  CI_prob <- c(pval/2, 0.5, 1-pval/2)
  iter <- dim(object$smp$h)[1]
  sel <- (burnin+1):iter
  tvec <- object$tvec
  
  # first, calculate summary statistics of h,H,S for each point
  # h
  h.q <- apply(object$smp$h[sel,,idx], c(2,3), quantile, probs=CI_prob)
  h.lower <- h.q[1,,]
  h.median <- h.q[2,,]
  h.upper <- h.q[3,,]
  # H
  H.q <- apply(object$smp$H[sel,,idx], c(2,3), quantile, probs=CI_prob)
  H.lower <- H.q[1,,]
  H.median <- H.q[2,,]
  H.upper <- H.q[3,,]
  # S
  S.q <- apply(object$smp$S[sel,,idx], c(2,3), quantile, probs=CI_prob)
  S.lower <- S.q[1,,]
  S.median <- S.q[2,,]
  S.upper <- S.q[3,,]
  
  if (popmean) {
    S.popmean <- apply(object$smp$S[,,idx], c(1,2), mean)
    S.popmean.q <- apply(S.popmean[sel,], 2, quantile, probs=CI_prob)
    S.popmean.lower <- S.popmean.q[1,]
    S.popmean.median <- S.popmean.q[2,]
    S.popmean.upper <- S.popmean.q[3,]
    
    if (make.plot) {
      S.popmean.ylim <- range(S.popmean.lower, S.popmean.median, S.popmean.upper)
      plot(tvec, S.popmean.median, type="l", xlab="Time", ylab="Population Survival Probability", ylim=S.popmean.ylim)
      lines(tvec, S.popmean.lower, lty=2)
      lines(tvec, S.popmean.upper, lty=2)
      lines(object$km.fit, col="red")
      legend("topright", legend = c("bsgw model", "kaplan-meyer"), col=c("black","red"), lty = c(1,1))
    }
  } else {
    S.popmean.lower <- NA
    S.popmean.median <- NA
    S.popmean.upper <- NA
  }
  
  # survival curves
  S.q <- apply(object$smp$S[sel,,idx], c(2,3), quantile, probs=CI_prob)
  S.lower <- S.q[1,,]
  S.median <- S.q[2,,]
  S.upper <- S.q[3,,]

  # pair-wise comparisons
  if (length(idx)==2) {
    if (tvec[1]==0) {
      tindex <- 2:length(tvec)
    } else {
      tindex <- 1:length(tvec)
    }
    idx1 <- idx[1]
    idx2 <- idx[2]
    # hazard ratio
    hr <- object$smp$h[,,idx2]/object$smp$h[,,idx1]
    hr.q <- apply(hr[sel,tindex], 2, quantile, probs=CI_prob)
    hr.lower <- hr.q[1,]
    hr.median <- hr.q[2,]
    hr.upper <- hr.q[3,]
    # survival diff
    S.diff <- object$smp$S[sel,,idx2]-object$smp$S[sel,,idx1]
    S.diff.q <- apply(S.diff, 2, quantile, probs=CI_prob)
    S.diff.lower <- S.diff.q[1,]
    S.diff.median <- S.diff.q[2,]
    S.diff.upper <- S.diff.q[3,]
    
    if (make.plot) {
      hr.range <- range(hr.lower, hr.median, hr.upper, 1.0)
      plot(tvec[tindex], hr.median, type="l", xlab="Time", ylab="Hazard Ratio"
           , ylim=hr.range, main=paste0("pval=, ", pval, ", idx1=", idx1, ", idx2=", idx2))
      lines(tvec[tindex], hr.lower, lty=2)
      lines(tvec[tindex], hr.upper, lty=2)
      abline(h=1.0, lty=2, col="red")
      
      S.range <- range(S.lower, S.median, S.upper)
      plot(tvec, S.median[,1], type="l", xlab="Time", ylab="Survival Probability", ylim=S.range, col="green"
           , main=paste0("pval=", pval))
      lines(tvec, S.lower[,1], lty=2, col="green")
      lines(tvec, S.upper[,1], lty=2, col="green")
      lines(tvec, S.median[,2], col="red")
      lines(tvec, S.lower[,2], lty=2, col="red")
      lines(tvec, S.upper[,2], lty=2, col="red")
      legend("topright", legend=c(paste0("idx1=",idx1),paste0("idx2=",idx2)), col=c("green","red"), lty=c(1,1))
      
      S.diff.range <- range(S.diff.lower, S.diff.median, S.diff.upper, 0.0)
      plot(tvec, S.diff.median, type="l", xlab="Time", ylab="Survival Probability Difference"
           , ylim=S.diff.range, main=paste0("pval=, ", pval, ", idx1=", idx1, ", idx2=", idx2))
      lines(tvec, S.diff.lower, lty=2)
      lines(tvec, S.diff.upper, lty=2)
      abline(h=0.0, lty=2, col="red")
    }
  } else {
    hr.lower <- NA
    hr.median <- NA
    hr.upper <- NA
    S.diff.lower <- NA
    S.diff.median <- NA
    S.diff.upper <- NA
  }
  
  return (list(lower=list(h=h.lower, H=H.lower, S=S.lower, hr=hr.lower, S.diff=S.diff.lower)
               , median=list(h=h.median, H=H.median, S=S.median, hr=hr.median, S.diff=S.diff.median)
               , upper=list(h=h.upper, H=H.upper, S=S.upper, hr=hr.upper, S.diff=S.diff.upper)
               , popmean=list(lower=list(S=S.popmean.lower), median=list(S=S.popmean.median), upper=list(S=S.popmean.upper))
               , km.fit=object$km.fit))
}

bsgw.crossval <- function(data, folds, all=FALSE, print.level=1, control=bsgw.control(), ncores=1, ...) {
  nfolds <- max(folds) # TODO: add more checks for validity of folds
  if (all) {
    if (ncores==1) {
      ret <- lapply (1:nfolds, function(n) {
        if (print.level>=1) cat("processing fold", n, "of", nfolds, "\n")
        flush.console()
        est <- bsgw(data=data[which(folds!=n),], control=control, print.level=print.level, ...)
        pred <- predict(est, newdata=data[which(folds==n),], burnin=control$burnin)
        ret <- max(pred$smp$loglike)
        attr(ret, "estobj") <- est
        return (ret)
      })
    } else { # parallel code; TODO: set upper bound on ncores based on maximum processor cores
      registerDoParallel(ncores)
      ret <- foreach (n=1:nfolds, .options.multicore=list(preschedule=FALSE)) %dopar% {
        if (print.level>=1) cat("processing fold", n, "of", nfolds, "\n")
        flush.console()
        est <- bsgw(data=data[which(folds!=n),], control=control, print.level=print.level, ...)
        pred <- predict(est, newdata=data[which(folds==n),], burnin=control$burnin)
        ret <- max(pred$smp$loglike)
        attr(ret, "estobj") <- est
        return (ret)
      }
    }
    fret <- sum(unlist(ret))
    estobjs <- list()
    for (n in 1:nfolds) estobjs[[n]] <- attr(ret[[n]], "estobj")
    attr(fret, "estobjs") <- estobjs
    return (fret)
  } else {
    if (ncores==1) {
      loglike <- sapply (1:nfolds, function(n) {
        if (print.level>=1) cat("processing fold", n, "of", nfolds, "\n")
        est <- bsgw(data=data[which(folds!=n),], control=control, print.level=print.level, ...)
        pred <- predict(est, newdata=data[which(folds==n),], burnin=control$burnin)
        return (max(pred$smp$loglike))
      })
      return (sum(loglike))
    } else {
      registerDoParallel(ncores)
      loglike <- foreach (n=1:nfolds, .options.multicore=list(preschedule=FALSE), .combine=c) %dopar% {
        if (print.level>=1) cat("processing fold", n, "of", nfolds, "\n")
        est <- bsgw(data=data[which(folds!=n),], control=control, print.level=print.level, ...)
        pred <- predict(est, newdata=data[which(folds==n),], burnin=control$burnin)
        return (max(pred$smp$loglike))
      }
      return (sum(loglike))
    }
  }
}

# TODO: under construction
bsgw.crossval.wrapper <- function(data, folds, all=FALSE, print.level=1, control=bsgw.control(), ncores=1
                                  , lambda.vec=exp(seq(from=log(0.01), to=log(100), length.out = 10))
                                  , lambdas.vec=NULL
                                  , lambda2=if (is.null(lambdas.vec)) cbind(lambda=lambda.vec, lambdas=lambda.vec)
                                  else as.matrix(expand.grid(lambda=lambda.vec, lambdas=lambdas.vec)), plot=TRUE, ...) {
  # TODO: need to impose upper bound on ncores by physical number of cores on system
  nfold <- length(unique(folds))
  # determining if to parallelize inner loop (over folds) or outer loop (over lambda's)
  # preference is inner due to better load balancing, unless we have more cores than folds
  # TODO: this logic can be improved, e.g. what if ncores=2 and nfold=3? then parallelization gain is only 3/2=1.5
  # but if there are many lambda combinations in outer loop then despite potentially uneven times, we can use dynamic
  # scheduling and gain better performance
  if (ncores<=nfold && nfold%%ncores==0) {
    if (print.level>=1) cat("using inner parallelization\n")
    ncores.inner <- ncores
    ncores.outer <- 1
  } else {
    if (print.level>=1) cat("using outer parallelization\n")
    ncores.inner <- 1
    ncores.outer <- ncores
  }
  nlambda <- nrow(lambda2)
  loglike <- rep(NA, nlambda)
  estobjs <- list()
  if (print.level>=1) cat("number of lambda combinations to test:", nlambda, "\n")
  if (ncores.outer==1) {
    for (i in 1:nlambda) {
      if (print.level>=1) cat("processing lambda combo", i, "of", nlambda, "\n")
      flush.console()
      control$lambda <- lambda2[i,"lambda"]
      control$lambdas <- lambda2[i,"lambdas"]
      ret <- bsgw.crossval(data=data, folds=folds, all=all, print.level=print.level, control=control, ncores=ncores.inner, ...)
      loglike[i] <- ret
      if (all) estobjs[[i]] <- attr(ret, "estobjs")
    }
  } else {
    registerDoParallel(ncores.outer)
    ret.all <- foreach (i=1:nlambda, .options.multicore=list(preschedule=FALSE)) %dopar% {
      if (print.level>=1) cat("processing lambda combo", i, "of", nlambda, "\n")
      flush.console()
      control$lambda <- lambda2[i,"lambda"]
      control$lambdas <- lambda2[i,"lambdas"]
      ret <- bsgw.crossval(data=data, folds=folds, all=all, print.level=print.level, control=control, ncores=ncores.inner, ...)
      return (list(loglike=ret, estobj=attr(ret, "estobjs")))
    }
    for (i in 1:nlambda) {
      loglike[i] <- ret.all[[i]]$loglike
      estobjs[[i]] <- ret.all[[i]]$estobj
    }
  }
  opt.index <- which(loglike==max(loglike))[1]
  lambda.opt <- lambda2[opt.index,"lambda"]
  lambdas.opt <- lambda2[opt.index,"lambdas"]
  
  ret <- list(lambda=lambda.opt, lambdas=lambdas.opt)
  attr(ret, "loglike.vec") <- loglike
  attr(ret, "loglike.opt") <- max(loglike)
  attr(ret, "lambda2") <- lambda2
  if (all) attr(ret, "estobjs") <- estobjs
  if (print.level>=1) cat("selected lambda's: ", c(lambda.opt, lambdas.opt), "\n")
  
  if (plot) {
    # we only plot when lambda and lambdas are the same, otherwise we need a 2D plot
    if (identical(attr(ret,"lambda2")[,"lambda"], attr(ret, "lambda2")[,"lambdas"])) {
      plot(lambda2[,"lambda"], loglike, type="l", log="x", xlab="Shrinkage", ylab="Log-likelihood")
    }
  }
  
  return (ret)
}





