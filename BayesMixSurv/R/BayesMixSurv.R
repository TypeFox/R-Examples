logLik.bayesmixsurv <- function(object, ...) {
  ret <- max(object$smp$loglike)
  # TODO: incomplete, how to find df in the presence of lambda? can we assume lambda only has an additive contribution
  # and therefore it can be ignored when comparing models of same lambda?
  class(ret) <- "logLik"
  return (ret)
}

# estimate
bayesmixsurv <- function(formula1, data, formula2=formula1, stratCol=NULL, weights, subset
  , na.action=na.fail, control=bayesmixsurv.control(), print.level=2) {
  # TODO: implement weights, subset, na.action
  mycall <- match.call()
  if (!missing(weights)) warning("weights argument not supported yet, this argument will be ignored")
  if (!missing(subset)) warning("subset argument not supported yet, this argument will be ignored")
  if (!identical(na.action,na.fail)) stop("na.action argument not supported yet; only na.fail is currently accepted")

  # TODO: need to make sure number of rows in X and X2 ends up being equal (watch out for removal of missing rows)
  mf1 <- model.frame(formula1, data, drop.unused.levels=TRUE, na.action = na.fail) # incorporate na.action argument
  mt1 <- attr(mf1, "terms")
  X1 <- model.matrix(mt1, mf1)
  y <- model.response(mf1, "numeric")
  colnamesX1 <- colnames(X1)
  if (colnamesX1[1]!="(Intercept)") stop("intercept term must be included in formula1")

  mf2 <- model.frame(formula2, data, drop.unused.levels=TRUE, na.action = na.fail) # incorporate na.action argument
  mt2 <- attr(mf2, "terms")
  X2 <- model.matrix(mt2, mf2)
  colnamesX2 <- colnames(X2)
  if (colnamesX2[1]!="(Intercept)") stop("intercept term must be included in formula2")

  if (is.list(control$scalex)) {
    X1 <- bayesmixsurv.scale(X1, apply.sc=control$scalex$apply.scale.X1, center=control$scalex$centerVec.X1, scale=control$scalex$scaleVec.X1)
    X2 <- bayesmixsurv.scale(X2, apply.sc=control$scalex$apply.scale.X2, center=control$scalex$centerVec.X2, scale=control$scalex$scaleVec.X2)
    control$scalex <- TRUE
  } else if (control$scalex) {
    X1 <- bayesmixsurv.scale(X1)
    X2 <- bayesmixsurv.scale(X2)
  }
  apply.scale.X1 <- attr(X1, "apply.scale")
  apply.scale.X2 <- attr(X2, "apply.scale")
  centerVec.X1 <- attr(X1, "centerVec")
  centerVec.X2 <- attr(X2, "centerVec")
  scaleVec.X1 <- attr(X1, "scaleVec")
  scaleVec.X2 <- attr(X2, "scaleVec")
  
  if (!is.null(stratCol)) {
    if ((stratCol %in% all.vars(formula1)) || (stratCol %in% all.vars(formula2))) stop("stratification column cannot be part of either formulas")
    data$stratColFactor <- as.factor(data[,stratCol]) # TODO: check that stratColFactor can be used as a name
    gformula <- ~ stratColFactor
    mfg <- model.frame(gformula, data, drop.unused.levels=TRUE)
    mtg <- attr(mfg, "terms")
    Xg <- model.matrix(mtg, mfg) # TODO: make sure intercept is first column
    stratContrasts <- attr(Xg, "contrasts")
    Xg <- Xg[,2:ncol(Xg),drop=F]
    stratXlevels <- .getXlevels(mtg, mfg)
    stratTerms <- mtg
    colnamesXg <- colnames(Xg)
    #Xg.debug <<- Xg
  } else {
    Xg <- NULL
    stratContrasts <- NULL
    stratXlevels <- NULL
    stratTerms <- NULL
    colnamesXg <- NULL
  }
  
  ret <- list(call=mycall, formula1=formula1, formula2=formula2, weights=rep(1,nrow(X1)), subset=1:nrow(X1)
              , na.action=na.action, control=control, X1=X1, X2=X2, y=y
              , contrasts1=attr(X1, "contrasts"), contrasts2=attr(X2, "contrasts")
              , xlevels1=.getXlevels(mt1, mf1), xlevels2=.getXlevels(mt2, mf2)
              , terms1=mt1, terms2=mt2
              , colnamesX1=colnamesX1, colnamesX2=colnamesX2
              , apply.scale.X1=apply.scale.X1, apply.scale.X2=apply.scale.X2
              , centerVec.X1=centerVec.X1, centerVec.X2=centerVec.X2
              , scaleVec.X1=scaleVec.X1, scaleVec.X2=scaleVec.X2
              , Xg=Xg, stratContrasts=stratContrasts, stratXlevels=stratXlevels, stratTerms=stratTerms, colnamesXg=colnamesXg
  )

  mcmc <- bayesmixsurv.mcmc(X1, X2, y[,1], y[,2], Xg, control$lambda1, control$lambda2, control$iter, control$single
                            , control$alpha2.fixed, control$alpha.boundary, control$sd.thresh, print.level, control$nskip)
  sel <- (control$burnin+1):control$iter
  median <- list(alpha1=median(mcmc$alpha1[sel]), alpha2=median(mcmc$alpha2[sel]), 
                 beta1=apply(mcmc$beta1[sel,,drop=F], 2, median), beta2=apply(mcmc$beta2[sel,,drop=F], 2, median)
                 , gamma=apply(mcmc$gamma[sel,,drop=F], 2, median), sigma.gamma=median(mcmc$sigma.gamma[sel]))
  
  km.fit <- survfit(bayesmixsurv.strip.formula(formula1), data)

  ret <- c(ret, list(idx1=mcmc$idx1, idx2=mcmc$idx2, median=median, max=list(loglike=max(mcmc$loglike))
                     , smp=list(alpha1=mcmc$alpha1, alpha2=mcmc$alpha2, beta1=mcmc$beta1, beta2=mcmc$beta2
                                , loglike=mcmc$loglike, gamma=mcmc$gamma, sigma.gamma=mcmc$sigma.gamma)
                     , km.fit=km.fit, tmax=max(y[,1])))
  class(ret) <- "bayesmixsurv"
  return (ret)
}

# print
print.bayesmixsurv <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("component-1 formula:\n")
  print(x$formula1)
  if (!x$control$single) {
    cat("component 2 formula:\n")
    print(x$formula2)
  } else {
    cat("(no component-2 formula in single-component mode)\n")
  }
  cat("component-1 shrinkage:", x$control$lambda1, "\n")
  if (!x$control$single) {
    cat("component-2 shrinkage:", x$control$lambda2, "\n")
  } else {
    cat("(no component-2 shrinkage in single-component mode)\n")
  }
  cat("MCMC iterations:", x$control$iter, "\n")
  cat("burn-in (median calc):", x$control$burnin, "\n")
  cat("Threshold on stdev of covariates:", x$control$sd.thresh, "\n")
  cat("Model matrix is scaled:", x$control$scalex, "\n")
  cat("component-1 coefficients:\n")
  beta1 <- x$median$beta1; names(beta1) <- x$colnamesX1
  print(beta1)
  if (!x$control$single) {
    cat("component-2 coefficients:\n")
    beta2 <- x$median$beta2; names(beta2) <- x$colnamesX2
    print(beta2)
  } else {
    cat("(no component-2 coefficients in single-component mode)\n")
  }
  cat("component-1 shape parameter:", x$median$alpha1, "\n")
  if (!x$control$single) {
    cat("component-2 shape parameter:", x$median$alpha2, "\n")
  } else {
    cat("(no component-2 shape parameter in single-component mode)\n")
  }
  cat("number of observations:", nrow(x$X1), "\n")
}

# plot
plot.bayesmixsurv <- function(x, pval=0.05, burnin=round(x$control$iter/2), nrow=2, ncol=3, ...) {
  iter <- x$control$iter
  sel <- (burnin+1):x$control$iter
  nsel <- length(sel)
  CI_prob <- c(pval/2, 0.5, 1-pval/2)
  nplot_per_page <- nrow*ncol
  # determine number of beta coefficients
  nbeta1 <- ncol(x$smp$beta1)
  nbeta2 <- ncol(x$smp$beta2)
  npage_beta1 <- ceiling(nbeta1/nplot_per_page)
  npage_beta2 <- ceiling(nbeta2/nplot_per_page)
  
  ## loglike and logpost
  par(mfrow=c(2,1))
  plot(x$smp$loglike, type="l", xlab="Iteration", ylab="Log-likelihood", main="Log-likelihood")
  plot(x$smp$loglike[sel], type="l", xlab="Iteration", ylab="Log-likelihood", main="Log-likelihood, Post-Burnin")
  
  ## traceplots
  # scale coefficients
  # beta1
  beta1_q <- apply(x$smp$beta1[sel,,drop=F], 2, quantile, probs=CI_prob)
  beta1_lower <- beta1_q[1,]
  beta1_median <- beta1_q[2,]
  beta1_upper <- beta1_q[3,]
  for (n in 1:npage_beta1) {
    par(mfrow=c(nrow,ncol))
    offset <- (n-1)*nplot_per_page
    for (i in 1:nplot_per_page) {
      if (offset+i<=nbeta1) {
        beta1_ylim <- range(x$smp$beta1[,offset+i], 0.0)
        plot(x$smp$beta1[,offset+i], type="l", xlab="Iteration", ylab="Sample Value", ylim=beta1_ylim
             , main = paste("beta1[", x$colnamesX1[offset+i], "]", sep=""))
        abline(h = 0)
        lines(sel, rep(beta1_lower[offset+i], nsel), lty=2, col="red")
        lines(sel, rep(beta1_median[offset+i], nsel), lty=2, col="red")
        lines(sel, rep(beta1_upper[offset+i], nsel), lty=2, col="red")
      }
    }
  }
  # beta2
  if (!x$control$single) {
    beta2_q <- apply(x$smp$beta2[sel,,drop=F], 2, quantile, probs=CI_prob)
    beta2_lower <- beta2_q[1,]
    beta2_median <- beta2_q[2,]
    beta2_upper <- beta2_q[3,]
    for (n in 1:npage_beta2) {
      par(mfrow=c(nrow,ncol))
      offset <- (n-1)*nplot_per_page
      for (i in 1:nplot_per_page) {
        if (offset+i<=nbeta2) {
          beta2_ylim <- range(x$smp$beta2[,offset+i], 0.0)
          plot(x$smp$beta2[,offset+i], type="l", xlab="Iteration", ylab="Sample Value", ylim=beta2_ylim
               , main = paste("beta2[", x$colnamesX2[offset+i], "]", sep=""))
          abline(h = 0)
          lines(sel, rep(beta2_lower[offset+i], nsel), lty=2, col="red")
          lines(sel, rep(beta2_median[offset+i], nsel), lty=2, col="red")
          lines(sel, rep(beta2_upper[offset+i], nsel), lty=2, col="red")
        }
      }
    }
  }
  # shape parameters
  par(mfrow=c(1,2))
  # alpha1
  alpha1_q <- quantile(x$smp$alpha1[sel], probs=CI_prob)
  alpha1_lower <- alpha1_q[1]
  alpha1_median <- alpha1_q[2]
  alpha1_upper <- alpha1_q[3]
  alpha1_ylim <- range(x$smp$alpha1)
  plot(x$smp$alpha1, type="l", xlab="Iteration", ylab="Sample Value", ylim=alpha1_ylim
       , main = paste("alpha1", sep=""))
  abline(h = 0)
  abline(h = alpha1_lower, lty=2, col="red")
  abline(h = alpha1_median, lty=2, col="red")
  abline(h = alpha1_upper, lty=2, col="red")
  if (!x$control$single) {
    # alpha2
    alpha2_q <- quantile(x$smp$alpha2[sel], probs=CI_prob)
    alpha2_lower <- alpha2_q[1]
    alpha2_median <- alpha2_q[2]
    alpha2_upper <- alpha2_q[3]
    alpha2_ylim <- range(x$smp$alpha2)
    plot(x$smp$alpha2, type="l", xlab="Iteration", ylab="Sample Value", ylim=alpha2_ylim
         , main = paste("alpha2", sep=""))
    abline(h = 0)
    abline(h = alpha2_lower, lty=2, col="red")
    abline(h = alpha2_median, lty=2, col="red")
    abline(h = alpha2_upper, lty=2, col="red")
  }
  
  ## autocorrelation plots
  # beta1
  for (n in 1:npage_beta1) {
    par(mfrow=c(nrow,ncol))
    offset <- (n-1)*nplot_per_page
    for (i in 1:nplot_per_page) {
      if (offset+i<=nbeta1) {
        if ((offset+i) %in% x$idx1) {
          acf(x$smp$beta1[sel,offset+i], main=paste("beta1[", x$colnamesX1[offset+i], "]", sep=""))
        } else {
          bayesmixsurv.empty.plot(main=paste("beta1[", x$colnamesX1[offset+i], "]", sep=""))
        }
      }
    }
  }
  if (!x$control$single) {
    # beta2
    for (n in 1:npage_beta2) {
      par(mfrow=c(nrow,ncol))
      offset <- (n-1)*nplot_per_page
      for (i in 1:nplot_per_page) {
        if (offset+i<=nbeta2) {
          if ((offset+i) %in% x$idx2) {
            acf(x$smp$beta2[sel,offset+i], main=paste("beta2[", x$colnamesX2[offset+i], "]", sep=""))
          } else {
            bayesmixsurv.empty.plot(main=paste("beta2[", x$colnamesX2[offset+i], "]", sep=""))
          }
        }
      }
    }
  }
    
  ## histograms
  # beta1
  for (n in 1:npage_beta1) {
    par(mfrow=c(nrow,ncol))
    offset <- (n-1)*nplot_per_page
    for (i in 1:nplot_per_page) {
      if (offset+i<=nbeta1) {
        if ((offset+i) %in% x$idx1) {
          hist(x$smp$beta1[sel,offset+i], xlab="Sample Value"
               , main=paste("beta1[", x$colnamesX1[offset+i], "]", sep=""))
          abline(v = 0)
          abline(v = beta1_median[offset+i], lty=2, col="red")
          abline(v = beta1_lower[offset+i], lty=3, col="red")
          abline(v = beta1_upper[offset+i], lty=3, col="red")
        } else {
          bayesmixsurv.empty.plot(main=paste("beta1[", x$colnamesX1[offset+i], "]", sep=""))
        }
      }
    }
  }
  if (!x$control$single) {
    # beta2
    for (n in 1:npage_beta2) {
      par(mfrow=c(nrow,ncol))
      offset <- (n-1)*nplot_per_page
      for (i in 1:nplot_per_page) {
        if (offset+i<=nbeta2) {
          if ((offset+i) %in% x$idx2) {
            hist(x$smp$beta2[sel,offset+i], xlab="Sample Value"
                 , main=paste("beta2[", x$colnamesX2[offset+i], "]", sep=""))
            abline(v = 0)
            abline(v = beta2_median[offset+i], lty=2, col="red")
            abline(v = beta2_lower[offset+i], lty=3, col="red")
            abline(v = beta2_upper[offset+i], lty=3, col="red")
          } else {
            bayesmixsurv.empty.plot(main=paste("beta2[", x$colnamesX2[offset+i], "]", sep=""))
          }
        }
      }
    }
  }
}

#summary
summary.bayesmixsurv <- function(object, pval=0.05, burnin=object$control$burnin, ...) {
  iter <- object$control$iter
  CI_prob <- c(pval/2, 0.5, 1-pval/2)
  sel <- (burnin+1):iter
  
  # alpha1, alpha2
  alpha1_q <- quantile(object$smp$alpha1[sel], prob=CI_prob)
  alpha1_lb <- alpha1_q[1]; alpha1_med <- alpha1_q[2]; alpha1_ub <- alpha1_q[3]
  alpha2_q <- quantile(object$smp$alpha2[sel], prob=CI_prob)
  alpha2_lb <- alpha2_q[1]; alpha2_med <- alpha2_q[2]; alpha2_ub <- alpha2_q[3]
  alpha1_pval <- bayesmixsurv.calc.pval(object$smp$alpha1[sel], ref=1.0)
  alpha2_pval <- bayesmixsurv.calc.pval(object$smp$alpha2[sel], ref=1.0)
  coefficients_alpha <- matrix(c(alpha1_med, alpha1_lb, alpha1_ub, alpha1_pval, alpha2_med, alpha2_lb, alpha2_ub, alpha2_pval)
                               , ncol=4, byrow = TRUE)
  dimnames(coefficients_alpha) <- list(c("alpha1","alpha2"), c("Estimate", "Lower Bound", "Upper Bound", "P-val"))
  
  # beta1
  beta1_q <- apply(object$smp$beta1[sel,,drop=F], 2, quantile, probs=CI_prob)
  beta1_lb <- beta1_q[1,]
  beta1_med <- beta1_q[2,]
  beta1_ub <- beta1_q[3,]
  beta1_pval <- apply(object$smp$beta1[sel,,drop=F], 2, bayesmixsurv.calc.pval, ref=0.0)
  coefficients_beta1 <- as.matrix(cbind(beta1_med, beta1_lb, beta1_ub, beta1_pval))
  dimnames(coefficients_beta1) <- list(object$colnamesX1, c("Estimate", "Lower Bound", "Upper Bound", "P-val"))
  # beta2
  beta2_q <- apply(object$smp$beta2[sel,,drop=F], 2, quantile, probs=CI_prob)
  beta2_lb <- beta2_q[1,]
  beta2_med <- beta2_q[2,]
  beta2_ub <- beta2_q[3,]
  beta2_pval <- apply(object$smp$beta2[sel,,drop=F], 2, bayesmixsurv.calc.pval, ref=0.0)
  coefficients_beta2 <- as.matrix(cbind(beta2_med, beta2_lb, beta2_ub, beta2_pval))
  dimnames(coefficients_beta2) <- list(object$colnamesX2, c("Estimate", "Lower Bound", "Upper Bound", "P-val"))
  
  # gamma
  if (!is.null(object$Xg)) {
    gamma_q <- apply(object$smp$gamma[sel,,drop=F], 2, quantile, probs=CI_prob)
    gamma_lb <- gamma_q[1,]
    gamma_med <- gamma_q[2,]
    gamma_ub <- gamma_q[3,]
    gamma_pval <- apply(object$smp$gamma[sel,,drop=F], 2, bayesmixsurv.calc.pval, ref=0.0)
    coefficients_gamma <- as.matrix(cbind(gamma_med, gamma_lb, gamma_ub, gamma_pval))
    dimnames(coefficients_gamma) <- list(object$colnamesXg, c("Estimate", "Lower Bound", "Upper Bound", "P-val"))
  } else {
    coefficients_gamma <- NULL
  }
  
  ret <- list(call=object$call, pval=pval, burnin=burnin, single=object$control$single
              , coefficients=list(alpha=coefficients_alpha, beta1=coefficients_beta1, beta2=coefficients_beta2
                                  , gamma=coefficients_gamma)
  )
  class(ret) <- "summary.bayesmixsurv"
  return (ret)
}

print.summary.bayesmixsurv <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("number of burn-in iterations discarded:", x$burnin, "\n")
  cat("confidence interval:", x$pval, "\n")
  cat("## shape coefficients ##\n")
  print(x$coefficients$alpha)
  cat("## scale coefficients ##\n")
  if (x$single) {
    print(x$coefficients$beta1)
  } else {
    cat("component 1:\n")
    print(x$coefficients$beta1)
    cat("component2:\n")
    print(x$coefficients$beta2)
  }
  if (!is.null(x$coefficients$gamma)) {
    cat("## stratification coefficients ##\n")
    print(x$coefficients$gamma)
  }
}

predict.bayesmixsurv <- function(object, newdata=NULL, tvec=NULL, burnin=object$control$burnin, ...) {
  iter <- object$control$iter
  alpha.min <- object$control$alpha.min
  alpha.max <- object$control$alpha.max
  
  tt1 <- object$terms1
  tt2 <- object$terms2
  Terms1 <- delete.response(tt1)
  Terms2 <- delete.response(tt2)
  
  if (is.null(newdata)) {
    nobs <- nrow(object$X1)
    X1 <- object$X1
    X2 <- object$X2
    Xg <- object$Xg
    km.fit <- object$km.fit
  } else {
    newdata <- droplevels(newdata)
    mf1 <- model.frame(Terms1, newdata, xlev = object$xlevels1)
    mf2 <- model.frame(Terms2, newdata, xlev = object$xlevels2)
    X1 <- model.matrix(Terms1, mf1, contrasts.arg = object$contrasts1)
    X2 <- model.matrix(Terms2, mf2, contrasts.arg = object$contrasts2)
    
    if (object$control$scalex) {
      X1 <- bayesmixsurv.scale(X1, apply.sc=object$apply.scale.X1, center=object$centerVec.X1, scale=object$scaleVec.X1)
      X2 <- bayesmixsurv.scale(X2, apply.sc=object$apply.scale.X2, center=object$centerVec.X2, scale=object$scaleVec.X2)
    }
    nobs <- nrow(newdata)
    
    if (!is.null(object$stratCol)) {
      newdata$stratColFactor <- as.factor(newdata[,object$stratCol])
      ttg <- object$stratTerms
      Termsg <- delete.response(ttg)
      
      index_seen_levels <- which(newdata$stratColFactor %in% object$stratXlevels$stratColFactor)
      index_unseen_levels <- setdiff(1:nrow(newdata), index_seen_levels)
      if (length(index_unseen_levels)>0) {
        mfg_seen <- model.frame(Termsg, newdata[index_seen_levels,], xlev = object$stratXlevels)
        Xg_seen <- model.matrix(Termsg, mfg_seen, contrasts.arg = object$stratContrasts)
        Xg <- matrix(NA, nrow=nrow(newdata), ncol=ncol(Xg_seen)); colnames(Xg) <- colnames(Xg_seen)
        Xg[index_seen_levels,] <- Xg_seen
        Xg[index_unseen_levels,] <- 1/ncol(Xg_seen)
      } else {
        mfg <- model.frame(Termsg, newdata, xlev = object$stratXlevels)
        Xg <- model.matrix(Termsg, mfg, contrasts.arg = object$stratContrasts)
      }
      Xg <- Xg[,-1] # dropping intercept term
    } else {
      Xg <- NULL
    }
  }
  
  if (!is.null(tvec)) {
    
    # TODO: we need an upper bound on length of tvec to avoid memory blow-up
    if (length(tvec)==1) tvec <- seq(from=0.0, to=object$tmax, length.out=tvec) # tvec is interpreted as number of time points
    
    nt <- length(tvec)
    tvec <- as.matrix(tvec)
    t_mat <- tvec[,rep(1,nobs)]
    
    ret <- lapply(1:iter, FUN=function(i) {
      xbeta1 <- X1%*%object$smp$beta1[i,]
      xbeta2 <- X2%*%object$smp$beta2[i,]
      alpha1 <- object$smp$alpha1[i]
      alpha2 <- object$smp$alpha2[i]
      
      exbeta1 <- as.matrix(exp(xbeta1))
      exbeta2 <- as.matrix(exp(xbeta2))
      exbeta1_mat <- t(exbeta1[,rep(1,nt)])
      exbeta2_mat <- t(exbeta2[,rep(1,nt)])
      
      H1tmp <- (t_mat^alpha1)*exbeta1_mat
      H2tmp <- (t_mat^alpha2)*exbeta2_mat
      
      h1tmp <- alpha1*(t_mat^(alpha1-1))*exbeta1_mat
      h2tmp <- alpha2*(t_mat^(alpha2-1))*exbeta2_mat
      
      return (list(h1=h1tmp, h2=h2tmp, H1=H1tmp, H2=H2tmp))
    })
    h1 <- array(NA, dim=c(iter, nt, nobs))
    h2 <- array(NA, dim=c(iter, nt, nobs))
    H1 <- array(NA, dim=c(iter, nt, nobs))
    H2 <- array(NA, dim=c(iter, nt, nobs))
    for (i in 1:iter) {
      h1[i,,] <- ret[[i]]$h1
      h2[i,,] <- ret[[i]]$h2
      H1[i,,] <- ret[[i]]$H1
      H2[i,,] <- ret[[i]]$H2
    }
    h <- h1+h2
    H <- H1+H2
    S <- exp(-H)
    
  } else {
    h1 <- NA
    h2 <- NA
    h <- NA
    H1 <- NA
    H2 <- NA
    H <- NA
    S <- NA
  }
  
  if (is.null(newdata)) {
    y <- object$y
    do_loglike <- T
  } else {
    Rterms <- drop.terms(tt1)
    if (all(all.vars(Rterms)[1:2] %in% colnames(newdata))) {
      mfy <- model.frame(Rterms, newdata, xlev = object$xlevels1) # TODO: add check to make sure response variable is available for newdata
      y <- model.response(mfy, "numeric")
      do_loglike <- T
      km.fit <- survfit(bayesmixsurv.strip.formula(object$formula1), newdata)
    } else {
      do_loglike <- F
      km.fit <- NULL
    }
  }
  
  if (do_loglike) { # TODO: include logpost
    loglike <- sapply(1:iter, FUN=function(i) {
      condprob.data(object$smp$alpha1[i], object$smp$beta1[i,], object$smp$alpha2[i], object$smp$beta2[i,], X1, X2, object$smp$gamma[i,], Xg, y[,1], y[,2])
    })
    loglike_median <- median(loglike[(burnin+1):iter])
  } else {
    loglike <- NA
    loglike_median <- NA
  }
  ret <- list(tvec=as.vector(tvec), burnin=burnin, median=list(loglike=loglike_median)
              , smp=list(h1=h1, h2=h2, h=h, H1=H1, H2=H2, H=H, S=S, loglike=loglike)
              , km.fit=km.fit
              )
  class(ret) <- "predict.bayesmixsurv"
  return (ret)
}

# summary of predict
summary.predict.bayesmixsurv <- function(object, idx=1:dim(object$smp$h)[3], burnin=object$burnin
                                 , pval=0.05, popmean=identical(idx,1:dim(object$smp$h)[3])
                                 , make.plot=TRUE, ...) {
  if (!all(idx %in% 1:dim(object$smp$h)[3])) {
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
      legend("topright", legend = c("bayesmixsurv model", "kaplan-meyer"), col=c("black","red"), lty = c(1,1))
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

# cross-validated log-likelihood
bayesmixsurv.crossval <- function(data, folds, all=FALSE, print.level=1, control=bayesmixsurv.control(), ...) {
  nfolds <- max(folds) # TODO: add more checks for validity of folds
  if (all) {
    ret <- lapply (1:nfolds, function(n) {
      if (print.level>=1) cat("processing fold", n, "of", nfolds, "\n")
      flush.console()
      est <- bayesmixsurv(data=data[which(folds!=n),], control=control, print.level=print.level, ...)
      pred <- predict(est, newdata=data[which(folds==n),], burnin=control$burnin)
      ret <- max(pred$smp$loglike)
      attr(ret, "estobj") <- est
      return (ret)
    })
    fret <- sum(unlist(ret))
    estobjs <- list()
    for (n in 1:nfolds) estobjs[[n]] <- attr(ret[[n]], "estobj")
    attr(fret, "estobjs") <- estobjs
    return (fret)
  } else {
    loglike <- sapply (1:nfolds, function(n) {
      if (print.level>=1) cat("processing fold", n, "of", nfolds, "\n")
      est <- bayesmixsurv(data=data[which(folds!=n),], control=control, print.level=print.level, ...)
      pred <- predict(est, newdata=data[which(folds==n),], burnin=control$burnin)
      return (max(pred$smp$loglike))
    })
    return (sum(loglike))
  }
}

bayesmixsurv.crossval.wrapper <- function(data, folds, all=FALSE, print.level=1, control=bayesmixsurv.control()
                                  , lambda.min=0.01, lambda.max=100, nlambda=10
                                  , lambda1.vec=exp(seq(from=log(lambda.min), to=log(lambda.max), length.out = nlambda))
                                  , lambda2.vec=NULL
                                  , lambda12=if (is.null(lambda2.vec)) cbind(lambda1=lambda1.vec, lambda2=lambda1.vec)
                                  else as.matrix(expand.grid(lambda1=lambda1.vec, lambda2=lambda2.vec)), plot=TRUE, ...) {
  nlambda <- nrow(lambda12)
  loglike <- rep(NA, nlambda)
  estobjs <- list()
  if (print.level>=1) cat("number of lambda combinations to test:", nlambda, "\n")
  for (i in 1:nlambda) { # TODO: keep this loop sequential due to load imbalance across different lambda's; parallelize folds within
    if (print.level>=1) cat("processing lambda combo", i, "of", nlambda, "\n")
    flush.console()
    control$lambda1 <- lambda12[i,"lambda1"]
    control$lambda2 <- lambda12[i,"lambda2"]
    ret <- bayesmixsurv.crossval(data=data, folds=folds, all=all, print.level=print.level, control=control, ...)
    loglike[i] <- ret
    if (all) estobjs[[i]] <- attr(ret, "estobjs")
  }
  opt.index <- which(loglike==max(loglike))[1]
  lambda1.opt <- lambda12[opt.index,"lambda1"]
  lambda2.opt <- lambda12[opt.index,"lambda2"]
  
  ret <- list(lambda1=lambda1.opt, lambda2=lambda2.opt)
  attr(ret, "loglike.vec") <- loglike
  attr(ret, "loglike.opt") <- max(loglike)
  attr(ret, "lambda12") <- lambda12
  if (all) attr(ret, "estobjs") <- estobjs
  if (print.level>=1) cat("selected lambda's: ", c(lambda1.opt, lambda2.opt), "\n")
  
  if (plot) {
    # we only plot when lambda and lambdas are the same, otherwise we need a 2D plot
    if (identical(attr(ret,"lambda12")[,"lambda1"], attr(ret, "lambda12")[,"lambda2"])) {
      plot(lambda12[,"lambda1"], loglike, type="l", log="x", xlab="Shrinkage", ylab="Log-likelihood")
    }
  }
  
  return (ret)
}

