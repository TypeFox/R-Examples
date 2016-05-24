#############################################################
#                                                           #
#	wle.glm function                                    #
#	Author: Claudio Agostinelli and Fatemah Alqallaf    #
#	E-mail: claudio@unive.it                            #
#	Date: March, 11, 2010                               #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                      and Fatemah Alqallaf                 #
#                                                           #
#############################################################

## Function developed from 'glm' R version 2.6.0
wle.glm <- function (formula, family = binomial, data, weights, subset, na.action, start = NULL, etastart, mustart, offset, control = list(glm=glm.control(...), wle=wle.glm.control()), model = TRUE, method = "wle.glm.fit", x = FALSE, y = TRUE, contrasts = NULL, dist.method="euclidean", ...) {
  call <- match.call()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  switch(method, model.frame = return(mf), wle.glm.fit = 1, stop("invalid 'method': ", method))
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(weights) && any(weights < 0)) 
    stop("negative weights not allowed")
  if (!is.null(offset)) {
    if (length(offset) == 1) 
      offset <- rep(offset, NROW(Y))
    else if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
  }
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")

### check the wle.glm.control entries
  if (is.null(control$wle$group))
    control$wle$group <- max(round(NROW(Y)/2), NCOL(X)+1)
  if (control$wle$group > NROW(Y))
    stop("'group' parameter must be not greater than NROW(Y)")
  maxboot <- sum(log(1:NROW(Y)))-(sum(log(1:control$wle$group))+sum(log(1:(NROW(Y)-control$wle$group))))
  if (log(control$wle$boot) > maxboot)
    stop("Bootstrap replication not in the range")
  if (is.null(control$wle$smooth) & !(family$family=='poisson' | family$family=='binomial')) {
    stop("'smooth parameter is missing without default")
  } 
  tot.sol <- 0
  not.conv <- 0
  iboot <- 0
  fit.store <- list()
  while (tot.sol < control$wle$num.sol & iboot < control$wle$boot) {
    iboot <- iboot + 1
    index.boot <- sample(1:NROW(Y), size=control$wle$group, replace=FALSE)
    if (NCOL(X) > 1)
      X.boot <- X[index.boot,]
    else
      X.boot <- X[index.boot]
    if (NCOL(Y) > 1)
      Y.boot <- Y[index.boot,]
    else
      Y.boot <- Y[index.boot]
  # initial values

    fit.boot <- try(glm.fit(x=X.boot, y=Y.boot, weights = weights[index.boot], start = start, etastart = etastart[index.boot], mustart = mustart[index.boot], offset = offset, family = family, control = control$glm, intercept = attr(mt, "intercept") > 0))
    if (is.list(fit.boot) && !any(is.na(fit.boot$coefficients))) {
      fit.boot$x <- data.frame(X.boot)
      fit.boot$y <- Y.boot
      fit.boot$terms <- mt
      fit.boot$family <- family
      class(fit.boot) <- "glm"
      fitted.Y <- try(predict.glm(fit.boot, newdata = if(is.null(data)) data.frame(X) else data, type = "response"))      
      dispersion <- if (family$family %in% c("poisson", "binomial")) {
        est.disp <- FALSE
        1
      } else if (fit.boot$df.residual > 0) {
        est.disp <- TRUE
        if (any(fit.boot$weights == 0)) 
          warning("observations with zero weight not used for calculating dispersion")
        if (control$wle$mle.dispersion) {
          if (family$family=="Gamma") { 
            wle.gamma.shape.glm(y=Y.boot, mu=fit.boot$fitted.values, deviance=fit.boot$deviance, df.residual=fit.boot$df.residual, prior.weights=weights[index.boot], wle.weights=rep(1, length(Y.boot)), it.lim=control$glm$maxit, eps.max=control$glm$epsilon, verbose=control$wle$verbose)$dispersion
          } else if (family$family=="inverse.gaussian") {
            wle.inversegaussian.lambda.glm(y=Y.boot, mu=fit.boot$fitted.values, prior.weights=weights[index.boot], wle.weights=rep(1, length(Y.boot)))$dispersion
          } else {
             stop("'WLE' method is not available for the estimation of the dispersion parameter in this family, set control$wle$mle.dispersion=FALSE to use the weighted method of moments")
          }
        } else {
          sum((fit.boot$weights * fit.boot$residuals^2)[fit.boot$weights > 0])/fit.boot$df.residual
        }
      } else {
        est.disp <- FALSE
        1
      }
###      cat('At the boot time ',dispersion, '\n')
  # initial wle.weights
      wle.weights <- wle.glm.weights(Y, X, fitted.values=fitted.Y, family = family, dispersion=dispersion, raf=control$wle$raf, tau=control$wle$tau, smooth=control$wle$smooth, asy.smooth=control$wle$asy.smooth, window.size=control$wle$window.size, use.asymptotic=control$wle$use.asymptotic, use.smooth=control$wle$use.smooth, tol=control$wle$tol, dist.method=dist.method, cutpoint=control$wle$cutpoint, powerdown=control$wle$powerdown)$weights
    
      fit <- wle.glm.fit(x = X, y = Y, weights = weights, wle.weights=wle.weights, start = fit.boot$coefficients, etastart = etastart, mustart = fit.boot$fitted.values, offset = offset, family = family, control = control, intercept = attr(mt, "intercept") > 0, dispersion=dispersion)
      if (length(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- wle.glm.fit(x = X[, "(Intercept)", drop = FALSE], y = Y, weights = weights, wle.weights=wle.weights, offset = offset, family = family, control = control, intercept = TRUE, dispersion=dispersion)$deviance
      }
      if (fit$converged) {
        fit$wle.weights <- fit$wle.weights/max(fit$wle.weights)
        if (tot.sol==0) {
          fit.store[[1]] <- fit
          w.store <- matrix(fit$wle.weights, nrow=1)
          tot.sol <- 1
        } else {
#####        cat(min(apply(abs(w.store-matrix(rep(fit$wle.weights, NROW(w.store)), ncol=NROW(Y), byrow=TRUE)), 1, max)), '\n')
          if (min(apply(abs(w.store-matrix(rep(fit$wle.weights, NROW(w.store)), ncol=NROW(Y), byrow=TRUE)), 1, max))>control$wle$equal) {
            tot.sol <- tot.sol+1
            fit.store[[tot.sol]] <- fit
            w.store <- rbind(w.store,fit$wle.weights)
          }
        }
      } else not.conv <- not.conv + 1
    } else {
      not.conv <- not.conv + 1
    }
  }
  if (tot.sol==0) {
    fit.store[[1]] <- fit
  }
    names(fit.store) <- paste('root', 1:length(fit.store), sep='')
  if (model) 
    fit.store$model <- mf
  fit.store$na.action <- attr(mf, "na.action")
  if (x)
    fit.store$x <- X
  if (!y)
    fit.store$y <- NULL
  else
    fit.store$y <- Y
  fit <- c(fit.store, list(family=fit.store$root1$family, call = call, formula = formula, terms = mt, data = data, offset = offset, control = control, method = method, contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)), tot.sol=tot.sol, not.conv=not.conv)
  class(fit) <- "wle.glm"
  return(fit)
}

#############################################################
#                                                           #
#	wle.glm.fit function                                #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 15, 2009                            #
#	Version: 0.2-1                                      #
#                                                           #
#	Copyright (C) 2009 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.glm.fit <- function(x, y, weights = NULL, wle.weights = rep(1, NROW(y)), start = NULL, etastart = NULL,  mustart = NULL, offset = rep(0, NROW(y)), family = gaussian(), control = list(glm=glm.control(), wle=wle.glm.control()), dist.method='euclidean', intercept = TRUE, dispersion=NULL) {
  if (is.null(weights))
    weights <- rep.int(1, NROW(y))

  if (family$family=='binomial' & NCOL(y) == 1 & !all((y==0) | (y==1))) {
    if (any((weights-round(weights))!=0)) {
      stop("Please submit the response variable as two columns (first column successes, second column unsuccesses), use 'weights' only for weighting the observations not for the total")
    } else {
      y <- cbind(round(y*weights), weights-round(y*weights))      
      weights <- rep.int(1, NROW(y))
    }
  }
  if (!is.null(control$wle$use.asymptotic) && NCOL(y) > 1)
    stop("The use of asymptotic residuals for a two columns response variable is not implemented yet.")

  if (NCOL(y) > 1) {
    ny <- NROW(y)
    wle.weights <- as.vector(t(wle.weights))
    weights <- c(weights, weights)
    y01 <- rbind(cbind(y[,1], rep(0, ny)), cbind(rep(0, ny), y[,2])) 
    x01 <- rbind(x, x)    
  } else {
    if (is.null(control$wle$delta)) {
      y01 <- y
    } else {
      delta <- control$wle$delta
      pibar <- mean(y)
      pihat <-  max(delta, min(1-delta, pibar)) 
      delta0 <- pihat*delta/(1+delta)
      delta1 <- (1+pihat*delta)/(1+delta)
      y01 <- delta0+(delta1-delta0)*y
      ## if (verbose)
      ##   cat(y01, '\n')
    }
    x01 <- x
  }
  
  pesia <- wle.weights 
  pesib <- wle.weights+1+control$wle$tol
  iter <- 0
####    iterations over wle weights until convergence of the wle weights  
  while(any(abs(pesib-pesia) > control$wle$tol) & iter < control$wle$max.iter) {
    iter <- iter + 1
    pesib <- pesia
    res.glm.fit <- try(glm.fit(x=x01, y=y01, weights = pesia*weights, start = start, etastart = etastart,  mustart = mustart, offset = offset, family = family, control = control$glm, intercept = intercept))
    if (is.list(res.glm.fit)) {###FIX ME!!!!!
      nvars <- NCOL(x01)
      EMPTY <- nvars == 0
      n.ok <- sum(pesia) - sum(res.glm.fit$weights==0)
      rank <- if(EMPTY) 0 else res.glm.fit$qr$rank
      res.glm.fit$df.null <- n.ok - as.integer(intercept)
      res.glm.fit$df.residual <- n.ok - rank
###    cat('df.residual ', res.glm.fit$df.residual, '\n')
      dispersion <- if (family$family %in% c("poisson", "binomial")) {
         est.disp <- FALSE
         1
      } else if (res.glm.fit$df.residual > 0) {
        est.disp <- TRUE
        if (any(res.glm.fit$weights == 0)) 
          warning("observations with zero weight not used for calculating dispersion")
        if (control$wle$mle.dispersion) {
          if (family$family=="Gamma") {
            wle.gamma.shape.glm(y=y01, mu=res.glm.fit$fitted.values, deviance=res.glm.fit$deviance, df.residual=res.glm.fit$df.residual, prior.weights=weights, wle.weights=pesia, it.lim=control$glm$maxit, eps.max=control$glm$epsilon, verbose=control$wle$verbose, dispersion=dispersion)$dispersion
          } else if (family$family=="inverse.gaussian") {
            wle.inversegaussian.lambda.glm(y=y01, mu=res.glm.fit$fitted.values, prior.weights=weights, wle.weights=pesia)$dispersion
          } else {
            stop("'WLE' method is not available for the estimation of the dispersion parameter in this family, set control$wle$mle.dispersion=FALSE to use the weighted method of moments")
          }
        } else {
            sum((pesia * res.glm.fit$weights * res.glm.fit$residuals^2)[res.glm.fit$weights > 0])/res.glm.fit$df.residual
        }
      } else {
        est.disp <- FALSE
        1
      }
###    cat('In the iterations ',dispersion, '\n')
      res.weights <- wle.glm.weights(y=y, x=x, fitted.values=res.glm.fit$fitted.values, family=family, dispersion=dispersion, raf=control$wle$raf, tau=control$wle$tau, smooth=control$wle$smooth, asy.smooth=control$wle$asy.smooth, window.size=control$wle$window.size, use.asymptotic=control$wle$use.asymptotic, use.smooth=control$wle$use.smooth, tol=control$wle$tol, dist.method=dist.method, cutpoint=control$wle$cutpoint, powerdown=control$wle$powerdown)
####    plot(pesia, ylim=c(0,1))
      pesia <- res.weights$weights
    } else {
      iter <- control$wle$max.iter
      start <- NULL
    }
  }

  ## calculate df
  nvars <- NCOL(x)
  EMPTY <- nvars == 0
  nobs <- NROW(y)
  n.ok <- sum(pesia) - sum(res.glm.fit$weights==0)
  rank <- if(EMPTY) 0 else res.glm.fit$qr$rank
  res.glm.fit$df.null <- n.ok - as.integer(intercept)
  res.glm.fit$df.residual <- n.ok - rank
#  ## calculate AIC
#  aic <- family$aic
#  ## 
#  eval(family$initialize)
#  ### controllare passaggio dei parametri
#  res.glm.fit$aic <- aic(y, n, res.glm.fit$fitted.values, res.glm.fit$weights, res.glm.fit$deviance) + 2*rank
  
  if (iter == control$wle$max.iter) res.glm.fit$converged <- FALSE
  res.glm.fit$prior.weights <- weights
  res.glm.fit$wle.weights <- pesia
  res.glm.fit$wle.asymptotic <- res.weights$asymptotic
  res.glm.fit$dispersion <- dispersion
  return(res.glm.fit)
}

#############################################################
#                                                           #
#	wle.glm.weights function                            #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: May, 04, 2011                                 #
#	Version: 0.6-1                                      #
#                                                           #
#	Copyright (C) 2011 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.glm.weights <- function(y, x, fitted.values, family = gaussian(), dispersion=1, raf='GKL', tau=0.1, smooth=NULL, asy.smooth=0.031, window.size=NULL, use.asymptotic=NULL, use.smooth=TRUE, tol=10^(-6), dist.method='euclidean', cutpoint=0, powerdown=1) {

  ## internal function for the evaluation of the weights given the Pearson residuals
  pesi <- function(x, raf, tau=0.1) {
    gkl <- function(x, tau) {
      if (tau!=0)
        x <- log(tau*x+1)/tau
      return(x)
    }
    pwd <- function(x, tau) {
      if(tau==Inf)
        x <- log(x+1)
      else
        x <- tau*((x + 1)^(1/tau) - 1)
      return(x)
    }
    x <- ifelse(x < 1e-10, 0, x)
    ww <- switch(raf,
       ##park+basu+2003.pdf
       GKL = gkl(x, tau),
       ##lindsay+1994.pdf                 
       PWD = pwd(x, tau),
       HD =  2*(sqrt(x + 1) - 1) ,
       NED =  2 - (2 + x)*exp(-x) ,
       SCHI2 =  1-(x^2/(x^2 +2)) )

    if (raf!='SCHI2') {
      ww <- (ww + 1)/(x + 1)
    }
    ww[ww > 1] <- 1
    ww[ww < 0] <- 0
    ww[is.infinite(x)] <- 0
    return(ww)
  }

  sottopesa <- function(x, cutpoint=0, powerdown=1)
    ifelse(x > cutpoint | x==0, x, x^((cutpoint/x)^powerdown))

  ## internal function to check duplicate in the 'x' matrix
  myunique <- function(x, window.size=NULL) {
  ## x: a distance matrix
    nx <- nrow(x)
    index0 <- index1 <- list()
    nodup <- rep(TRUE, nx)
    j <- 0
    while (any(nodup)) {
      j <- j + 1
      jj <- min((1:nx)[nodup])
      ind0 <- (1:nx)[x[jj,] <= .Machine$double.eps*2]
      if (!is.null(window.size))
        ind1 <- (1:nx)[x[jj,] <= window.size]
      else
        ind1 <- ind0
      nodup[ind0] <- FALSE
      index0[[j]] <- ind0
      index1[[j]] <- ind1
    }
    result <- list(index0, index1)
    return(result)
  }

########## BINOMIAL AND POISSON
  if (family$family=='binomial' | family$family=='poisson') {
    ny <- NROW(y)
    if (NCOL(y)==1)
      weights <- rep(0, ny)
    else
      weights <- matrix(0, nrow=ny, ncol=2)
    asy <- rep(FALSE, ny)
    dx <- as.matrix(dist(x, method=dist.method))
    ressplit <- myunique(dx, window.size=window.size)
    splitx <- ressplit[[1]]
    splitxx <- ressplit[[2]]
    nny <- sapply(splitx, length)
    nnny <- sapply(splitxx, length)
    if (!is.null(use.asymptotic) && any(nnny <= use.asymptotic)) {
      residall <- residualsAnscombe(y=y, mu=fitted.values, family=family)
##        cat(summary(residall), '\n')
##        cat(sd(residall), '\n')
##        cat(mad(residall), '\n')
      ff.den <- density(residall, bw=asy.smooth, na.rm=TRUE)
##    plot(ff.den)
      ff.smoothed <- approxfun(x=ff.den$x, y=ff.den$y, rule=2)
      ffall <- ff.smoothed(v=residall)
      mmall <- dnorm(residall, mean=0, sd=sqrt(1+asy.smooth))
      ddall <- ffall/mmall - 1    
      asy.weights <- pesi(x=ddall, raf=raf, tau=tau)

##      cat('maxww', maxww, '\n')
    }

    for (j in 1:length(splitx)) {
      if (is.null(use.asymptotic) || nnny[j] > use.asymptotic) {
        asymptotic <- FALSE
        if (family$family=='binomial') {
          if (NCOL(y) == 1) {
            ytemp <- y[splitx[[j]]]
            ff <- rep(1, nny[j])
            tff <- table(y[splitxx[[j]]])
            tff <- tff/sum(tff)
            if (length(tff)==2) {
              ff[ytemp==0] <- tff[1]
              ff[ytemp==1] <- tff[2]
            }
            mm <- dbinom(ytemp,size=1,prob=fitted.values[splitx[[j]]][1])
          } else {
            ytemp <- y[splitx[[j]], , drop=FALSE]
            tff <- apply(y[splitxx[[j]], , drop = FALSE], 2, sum)
            tff <- tff/sum(tff)
            ff <- matrix(tff, nrow=nny[j], ncol=2, byrow=TRUE)
            mm <- dbinom(c(1,0), size=1, prob=fitted.values[splitx[[j]]][1])
            mm <- matrix(mm, nrow=nny[j], ncol=2, byrow=TRUE)
          }
        } else if (family$family=='poisson') {
          ytemp <- y[splitx[[j]]]
          ff <- rep(0,nny[j])
          tff <- table(y[splitxx[[j]]])
          tff <- tff/sum(tff)
          nff <- as.numeric(names(tff))
          for (i in 1:nny[j]) {
            ff[i] <- tff[nff==ytemp[i]] 
          }
          mm <- dpois(ytemp,lambda=fitted.values[splitx[[j]]][1])
        }
#        else {
#          stop('family not implemented yet')
#        } 
      } else {
        ##### quasi--family
        asymptotic <- TRUE
        rtemp <- residall[splitx[[j]]]
        ff <- ff.smoothed(v=rtemp)
        mm <- dnorm(rtemp, mean=0, sd=sqrt(1+asy.smooth))
      }
      dd <- ff/mm - 1

#    cat('ff', ff, '\n')
#    cat('mm', mm, '\n')
#    cat('dd', dd, '\n')
#    cat('ww', ww, '\n')
      if (NCOL(y)==1) {
        ww <- pesi(x=dd, raf=raf, tau=tau)
        weights[splitx[[j]]] <- ww
        asy[splitx[[j]]] <- asymptotic
      } else {
        ww <- apply(dd, 2, function(x) pesi(x, raf=raf, tau=tau))
        weights[splitx[[j]],] <- ww
      }
    }
    if (sum(asy) < ny & sum(asy) >= ny/2) {
      weights[asy] <- quantile((weights[!asy]/asy.weights[!asy]), 0.75)*weights[asy]
    } else if (sum(asy) > 0 & sum(asy) < ny/2) {
      weights[!asy] <- quantile((weights[asy]/asy.weights[asy]), 0.75)*weights[!asy]
    }
    weights[weights > 1] <- 1
########## QUASIBINOMIAL AND QUASIPOISSON    
  } else if (family$family=='quasibinomial' | family$family=='quasipoisson') {
    resid <- (y - fitted.values)/sqrt(family$variance(fitted.values))
    ff.den <- density(resid, bw=asy.smooth, na.rm=TRUE)
    ff.smoothed <- approxfun(x=ff.den$x, y=ff.den$y, rule=2)
    ffall <- ff.smoothed(v=resid)
    mmall <- dnorm(resid, mean=0, sd=sqrt(1+asy.smooth))
    ddall <- ffall/mmall - 1
    weights <- pesi(x=ddall, raf=raf, tau=tau)
    asy <- rep(TRUE, length(weights))
########## GAMMA ####################################    
  } else if (family$family=='Gamma') {
    ny <- length(y)
    weights <- rep(0, ny)
    asy <- rep(FALSE, ny)
    dx <- as.matrix(dist(x, method=dist.method))
    ressplit <- myunique(dx, window.size=window.size)
    splitx <- ressplit[[1]]
    splitxx <- ressplit[[2]]
    nny <- sapply(splitx, length)
    nnny <- sapply(splitxx, length)
    if (!is.null(use.asymptotic) && any(nnny <= use.asymptotic)) {
      residall <- residualsAnscombe(y=y, mu=fitted.values, family=family)/sqrt(dispersion)
      residall <- residall - median(residall)
##        cat(summary(residall), '\n')
##        cat(sd(residall), '\n')
##        cat(mad(residall), '\n')
      ff.den <- density(residall, bw=asy.smooth, na.rm=TRUE)
##    plot(ff.den)
      ff.smoothed <- approxfun(x=ff.den$x, y=ff.den$y, rule=2)
      ffall <- ff.smoothed(v=residall)
      mmall <- dnorm(residall, mean=0, sd=sqrt(1+asy.smooth))
      ddall <- ffall/mmall - 1    
      asy.weights <- pesi(x=ddall, raf=raf, tau=tau)
##      asy.weights[is.na(asy.weights)] <- 0
##      cat('maxww', maxww, '\n')
    }

    for (j in 1:length(splitx)) {
      if (is.null(use.asymptotic) || nnny[j] > use.asymptotic) {
        asymptotic <- FALSE
        ytemp <- y[splitx[[j]]]
        xtemp <- y[splitxx[[j]]]
        ### phi dispersion parameter
        ### mu mean parameter
        shape <- 1/dispersion ### 1/phi
        rate <- 1/(fitted.values[splitx[[j]]][1]*dispersion) ### 1/(mu*phi)
        temp <- (fitted.values[splitx[[j]]][1])^2*dispersion ### temp <- shape/rate^2

###        cat('Shape ', shape, '\n')
###        cat('Rate ', rate, '\n')
###        cat('temp ', temp, '\n')
        dsup <- max(xtemp)+ 3*smooth*temp
###        cat('dsup ', dsup, '\n')       
        tol.int <- tol*10^(-4)
        z <- .Fortran("wlegamma",
	    as.double(xtemp),
            as.double(ytemp),
	    as.integer(nnny[j]),
	    as.integer(nny[j]),                     
	    as.integer(1), ### Here the RAF is fixed to Hellinger,
                           ### hence do not use the weights from this function
                           ### but use the function pesi from the Pearson
                           ### Residuals
            as.double(1.0),
	    as.double(smooth*temp),
            as.integer(1*use.smooth),
            as.double(dsup),
	    as.double(tol),
            as.double(tol.int),
	    as.double(rate),
	    as.double(shape),
	    ww=double(nny[j]),
	    ff=double(nny[j]),
	    mm=double(nny[j]),
            PACKAGE = "wle")
        ff <- z$ff
        mm <- z$mm
      } else {
        ##### quasi--family
        asymptotic <- TRUE
        rtemp <- residall[splitx[[j]]]
        ff <- ff.smoothed(v=rtemp)
        mm <- dnorm(rtemp, mean=0, sd=sqrt(1+asy.smooth))
      }
      dd <- ff/mm - 1

#    cat('ff', ff, '\n')
#    cat('mm', mm, '\n')
#    cat('dd', dd, '\n')
      ww <- pesi(x=dd, raf=raf, tau=tau)
#    cat('ww', ww, '\n')

      weights[splitx[[j]]] <- ww
      asy[splitx[[j]]] <- asymptotic
    }
    if (sum(asy) < ny & sum(asy) >= ny/2) {
      weights[asy] <- quantile((weights[!asy]/asy.weights[!asy]), 0.75)*weights[asy]
    } else if (sum(asy) > 0 & sum(asy) < ny/2) {
      weights[!asy] <- quantile((weights[asy]/asy.weights[asy]), 0.75)*weights[!asy]
    }
    weights[is.na(weights)] <- 0
    weights[weights > 1] <- 1
    weights <- sottopesa(weights, cutpoint=cutpoint, powerdown=powerdown)
########## INVERSE GAUSSIAN ####################################    
  } else if (family$family=='inverse.gaussian') {
    ny <- length(y)
    weights <- rep(0, ny)
    asy <- rep(FALSE, ny)
    dx <- as.matrix(dist(x, method=dist.method))
    ressplit <- myunique(dx, window.size=window.size)
    splitx <- ressplit[[1]]
    splitxx <- ressplit[[2]]
    nny <- sapply(splitx, length)
    nnny <- sapply(splitxx, length)
    if (!is.null(use.asymptotic) && any(nnny <= use.asymptotic)) {
      residall <- residualsAnscombe(y=y, mu=fitted.values, family=family)
      residall <- residall - median(residall)
##        cat(summary(residall), '\n')
##        cat(sd(residall), '\n')
##        cat(mad(residall), '\n')
      ff.den <- density(residall, bw=asy.smooth, na.rm=TRUE)
##    plot(ff.den)
      ff.smoothed <- approxfun(x=ff.den$x, y=ff.den$y, rule=2)
      ffall <- ff.smoothed(v=residall)
      mmall <- dnorm(residall, mean=0, sd=sqrt(1+asy.smooth))
      ddall <- ffall/mmall - 1
      asy.weights <- pesi(x=ddall, raf=raf, tau=tau)
##      asy.weights[is.na(asy.weights)] <- 0
##      cat('maxww', maxww, '\n')
    }

    for (j in 1:length(splitx)) {
      if (is.null(use.asymptotic) || nnny[j] > use.asymptotic) {
        asymptotic <- FALSE
        ytemp <- y[splitx[[j]]]
        xtemp <- y[splitxx[[j]]]
        ### mu mean parameter
        mu <- fitted.values[splitx[[j]]][1]
###        cat('Mu ', mu, '\n')
        ### dispersion is 1/lambda
        dsup <- max(xtemp)+ 3*smooth*dispersion
###        cat('dsup ', dsup, '\n')       
        tol.int <- tol*10^(-4)
        z <- .Fortran("wleinvga",
	    as.double(xtemp),
            as.double(ytemp),
	    as.integer(nnny[j]),
	    as.integer(nny[j]),                     
	    as.integer(1), ### Here the RAF is fixed to Hellinger,
                           ### hence do not use the weights from this function
                           ### but use the function pesi from the Pearson
                           ### Residuals
            as.double(1.0),
	    as.double(smooth*dispersion),
            as.integer(1*use.smooth),
            as.double(dsup),
	    as.double(tol),
            as.double(tol.int),
	    as.double(mu),
	    as.double(1/dispersion),
	    ww=double(nny[j]),
	    ff=double(nny[j]),
	    mm=double(nny[j]),
            PACKAGE = "wle")
        ff <- z$ff
        mm <- z$mm
      } else {
        ##### quasi--family
        asymptotic <- TRUE
        rtemp <- residall[splitx[[j]]]
        ff <- ff.smoothed(v=rtemp)
        mm <- dnorm(rtemp, mean=0, sd=sqrt(1+asy.smooth))
      }
      dd <- ff/mm - 1

#    cat('ff', ff, '\n')
#    cat('mm', mm, '\n')
#    cat('dd', dd, '\n')
      ww <- pesi(x=dd, raf=raf, tau=tau)
#    cat('ww', ww, '\n')

      weights[splitx[[j]]] <- ww
      asy[splitx[[j]]] <- asymptotic
    }
    if (sum(asy) < ny & sum(asy) >= ny/2) {
      weights[asy] <- quantile((weights[!asy]/asy.weights[!asy]), 0.75)*weights[asy]
    } else if (sum(asy) > 0 & sum(asy) < ny/2) {
      weights[!asy] <- quantile((weights[asy]/asy.weights[asy]), 0.75)*weights[!asy]
    }
    weights[is.na(weights)] <- 0
    weights[weights > 1] <- 1
    weights <- sottopesa(weights, cutpoint=cutpoint, powerdown=powerdown)
    weights[weights < 0.01] <- 0
#######################################################################
  } else 
    stop("Only family 'gaussian', 'Gamma', 'inverse.gaussian', 'binomial', 'poisson', 'quasibinomial' and 'quasipoisson' are implemented (for now)")
  ###cat('WLE weights ', weights, '\n')
  result <- list(weights=weights, asymptotic=asy, grouped=nny, grouped.windows=nny)
  return(result)
}   

#############################################################
#                                                           #
#	print.wle.glm function                              #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: February, 25, 2010                            #
#	Version: 0.1-2                                      #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                                                           #
#############################################################
## from print.glm in R 2.6.0

print.wle.glm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (x$tot.sol > 0) {
      for (i in 1:x$tot.sol) {
        cat("Root: ", i, "\n")
        y <- extractRoot(x, root=i)
        print.wle.glm.root(x=y, digits=digits, ...)
      }
    } else {
      cat('No converged solutions were found\n\n')
    }
    cat("\n")
    cat("\nNumber of solutions ", x$tot.sol, "\n")
    cat("\n")    
    invisible(x)
}

###############################################################
print.wle.glm.root <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Coefficients")
  if (is.character(co <- x$contrasts)) 
    cat("  [contrasts: ", apply(cbind(names(co), co), 1, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits), 
    print.gap = 2, quote = FALSE)
  cat("\nThe following statistics are valid only if this model is the FULL model in a model selection procedure\n")
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", 
    x$df.residual, "Residual\n")
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
    cat("Null Deviance:\t   ", format(signif(x$null.deviance, 
      digits)), "\nResidual Deviance:", format(signif(x$deviance, 
      digits)), "\tAIC:", format(signif(x$aic, digits)), "\n")
    cat("Converged:\t   ", x$converged,"\n\n")
  invisible(x)
}

###############################################################
deviance.wle.glm <- function(object, root="all", ...) {
  if (root=="all")
    root <- 1:object$tot.sol
  root <- round(as.numeric(root))
  root <- root[root <= object$tot.sol & root > 0]
  if (!length(root))
    stop('not valid root numbers provided')
  dev <- rep(0, length(root))
  for (iroot in 1:length(root)) {
    dev[iroot] <- eval(parse(text=paste('object$root', root[iroot], '$deviance', sep='')))
  }
  names(dev) <- paste('root', root, sep='')
  return(dev)
}

## from effects.glm
effects.wle.glm <- function(object, root="all", ...) {
  if (root=="all")
    root <- 1:object$tot.sol
  root <- round(as.numeric(root))
  root <- root[root <= object$tot.sol & root > 0]
  if (length(root)==1) {
    eff <- eval(parse(text=paste('object$root', root, '$effects', sep='')))
  } else {
    eff <- vector(length=0)
    for (iroot in root) {
      eff <- rbind(eff, eval(parse(text=paste('object$root', iroot, '$effects', sep=''))))
    }
    rownames(eff) <- paste('root', root, sep='')
  }
  return(eff)
}

## from family.glm
family.wle.glm <- function(object, ...) object$family

## from model.frame.glm
model.frame.wle.glm <- function (formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call
        fcall$method <- "model.frame"
        fcall[[1L]] <- as.name("wle.glm")
        fcall[names(nargs)] <- nargs
#       env <- environment(fcall$formula)  # always NULL
        env <- environment(formula$terms)
        if (is.null(env)) env <- parent.frame()
        eval(fcall, env)
    }
    else formula$model
}

fitted.wle.glm <- function(object, root="all", ...) {
  if (root=="all")
    root <- 1:object$tot.sol
  root <- round(as.numeric(root))
  root <- root[root <= object$tot.sol & root > 0]
  if (!length(root))
    stop('not valid root numbers provided')
  if (length(root)==1) {
    fitted.values <- eval(parse(text=paste('object$root', root, '$fitted.values', sep='')))
  } else {
    fitted.values <- vector(length=0)
    for (iroot in root) {
      fitted.values <- rbind(fitted.values, eval(parse(text=paste('object$root', iroot, '$fitted.values', sep=''))))
    }
    rownames(fitted.values) <- paste('root', root, sep='')
  }
  return(fitted.values)
}

coef.wle.glm <- function(object, root="all", ...) {
  if (root=="all")
    root <- 1:object$tot.sol
  root <- round(as.numeric(root))
  root <- root[root <= object$tot.sol & root > 0]
  if (!length(root))
    stop('not valid root numbers provided')
  if (length(root)==1) {
    coeff <- eval(parse(text=paste('object$root', root, '$coefficients', sep='')))
  } else {
    coeff <- vector(length=0)
    for (iroot in root) {
      coeff <- rbind(coeff, eval(parse(text=paste('object$root', iroot, '$coefficients', sep=''))))
    }
    rownames(coeff) <- paste('root', root, sep='')
  }
  return(coeff)
}

## from weights.glm
weights.wle.glm <- function(object, type = c("prior", "working", "wle"), root="all", ...) {
  type <- match.arg(type)
  type <- switch(type,
            prior = "prior.weights",
            working = "weights",
            wle = "wle.weights")  
  if (root=="all")
    root <- 1:object$tot.sol
  root <- round(as.numeric(root))
  root <- root[root <= object$tot.sol & root > 0]
  if (!length(root))
    stop('not valid root numbers provided')
  if (length(root)==1) {
    res <- eval(parse(text=paste('object$root', root, '$', type, sep='')))
  } else {
    res <- vector(length=0)
    for (iroot in root) {
      res <- rbind(res, eval(parse(text=paste('object$root', iroot, '$', type, sep=''))))
    }
    rownames(res) <- paste('root', root, sep='')
  }
    
  if(is.null(object$na.action)) res
  else naresid(object$na.action, res)
}

## from formula.glm
formula.wle.glm <- function(x, ...) {
    form <- x$formula
    if( !is.null(form) ) {
        form <- formula(x$terms) # has . expanded
        environment(form) <- environment(x$formula)
        form
    } else formula(x$terms)
}

## from residuals.glm
residuals.wle.glm <- function(object, type = c("deviance", "pearson", "working","response", "partial"), root="all", ...) {
  type <- match.arg(type)
  if (type == "partial") 
    .NotYetImplemented()
  res.internal <- function(type, y, dev.resids, variance, mu.eta, r, mu, wts, eta, df.res) {
    switch(type,
           deviance=,pearson=,response=
           if(is.null(y)) {
               y <-  mu + r * mu.eta(eta)
           })
    res <- switch(type,
                  deviance = if(df.res > 0) {
                      d.res <- sqrt(pmax(dev.resids(y, mu, wts),
 0))
                      ifelse(y > mu, d.res, -d.res)
                  } else rep.int(0, length(mu)),
                  pearson = (y-mu)*sqrt(wts)/sqrt(variance(mu)),
                  working = r,
                  response = y - mu,
                  partial = r
                  )
    if(!is.null(object$na.action))
        res <- naresid(object$na.action, res)
    if (type == "partial") ## need to avoid doing naresid() twice.
        .NotYetImplemented() ##res <- res+predict(object, type="terms")
    return(res)
  }

  y <- object$y
  dev.resids <- object$family$dev.resids
  variance <- object$family$variance
  mu.eta <- object$family$mu.eta

  if (root=="all")
    root <- 1:object$tot.sol
  root <- round(as.numeric(root))
  root <- root[root <= object$tot.sol & root > 0]
  if (!length(root))
    stop('not valid root numbers provided')
  if (length(root)==1) {
    r <- eval(parse(text=paste('object$root', root, '$residuals', sep='')))
    mu <- eval(parse(text=paste('object$root', root, '$fitted.values', sep='')))
    wts <- eval(parse(text=paste('object$root', root, '$prior.weights', sep='')))
    eta <- eval(parse(text=paste('object$root', root, '$linear.predictors', sep='')))
    df.res <- eval(parse(text=paste('object$root', root, '$df.residual', sep='')))
    res <- res.internal(type, y, dev.resids, variance, mu.eta, r, mu, wts, eta, df.res)
  } else {
    res <- vector(length=0)
    for (iroot in root) {
      r <- eval(parse(text=paste('object$root', iroot, '$residuals', sep='')))
      mu <- eval(parse(text=paste('object$root', iroot, '$fitted.values', sep='')))
      wts <- eval(parse(text=paste('object$root', iroot, '$prior.weights', sep='')))
      eta <- eval(parse(text=paste('object$root', iroot, '$linear.predictors', sep='')))
      df.res <- eval(parse(text=paste('object$root', iroot, '$df.residual', sep='')))
      res <- rbind(res, res.internal(type, y, dev.resids, variance, mu.eta, r, mu, wts, eta, df.res))
    }
    rownames(res) <- paste('root', root, sep='')
  }
  return(res)
}

#################################
summary.wle.glm <- function(object, root=1, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...) {
  if (object$tot.sol>0) {
    root <- round(as.numeric(root))
    if (length(root)!=1)
      stop("'root' arguments must be a vector of length 1")
    if (root > object$tot.sol | root <= 0)
      stop("'root' must be between 1 and object$tot.sol")
    object <- extractRoot.wle.glm(object, root=root)
    res <- summary.wle.glm.root(object=object, dispersion=dispersion, correlation=correlation, symbolic.cor=symbolic.cor, ...)
    return(res)
  } else {
    cat('No converged solutions were found\n')
  }
}

## from summary.glm
summary.wle.glm.root <- function(object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...) {
    est.disp <- FALSE
    df.r <- object$df.residual
    if(is.null(dispersion))	# calculate dispersion if needed
	dispersion <-
	    if(object$family$family %in% c("poisson", "binomial"))  1
	    else if(df.r > 0) {
                est.disp <- TRUE
		if(any(object$weights==0))
		    warning("observations with zero weight not used for calculating dispersion")
                if (object$family$family=="Gamma" | object$family$family=="inverse.gaussian")
                  object$dispersion
                else 
  		  sum((object$wle.weights*object$weights*object$residuals^2)[object$weights > 0])/ df.r
	    } else {
                est.disp <- TRUE
                NaN
            }

    ## calculate scaled and unscaled covariance matrix

    aliased <- is.na(coef(object))  # used in print method
    p <- object$rank
    if (p > 0) {
        p1 <- 1L:p
        Qr <- object$qr
        ## WATCHIT! doesn't this rely on pivoting not permuting 1L:p? -- that's quaranteed
        coef.p <- object$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1,p1,drop=FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p),names(coef.p))
        covmat <- dispersion*covmat.unscaled
        var.cf <- diag(covmat)

        ## calculate coef table

        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err

        dn <- c("Estimate", "Std. Error")
        if(!est.disp) { # known dispersion
            pvalue <- 2*pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p),
                                         c(dn, "z value","Pr(>|z|)"))
        } else if(df.r > 0) {
            pvalue <- 2*pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p),
                                         c(dn, "t value","Pr(>|t|)"))
        } else { # df.r == 0
            coef.table <- cbind(coef.p, NaN, NaN, NaN)
            dimnames(coef.table) <- list(names(coef.p),
                                         c(dn, "t value","Pr(>|t|)"))
        }
        df.f <- NCOL(Qr$qr)
    } else {
        coef.table <- matrix(, 0L, 4L)
        dimnames(coef.table) <-
            list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0L, 0L)
        df.f <- length(aliased)
    }
    ## return answer

    ## these need not all exist, e.g. na.action.
    keep <- match(c("call","terms","family","deviance", "aic",
		      "contrasts", "df.residual","null.deviance","df.null",
                      "iter", "na.action"), names(object), 0L)
    ans <- c(object[keep],
	     list(deviance.resid = residuals(object, type = "deviance"),
		  coefficients = coef.table,
                  aliased = aliased,
		  dispersion = dispersion,
		  df = c(object$rank, df.r, df.f),
		  cov.unscaled = covmat.unscaled,
		  cov.scaled = covmat))

    if(correlation && p > 0) {
	dd <- sqrt(diag(covmat.unscaled))
	ans$correlation <-
	    covmat.unscaled/outer(dd,dd)
	ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.wle.glm"
    return(ans)
}

# from print.summary.glm
print.summary.wle.glm <- function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, signif.stars = getOption("show.signif.stars"), ...) {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")
    cat("Deviance Residuals: \n")
    if(x$df.residual > 5) {
	x$deviance.resid <- quantile(x$deviance.resid,na.rm=TRUE)
	names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
    }
    print.default(x$deviance.resid, digits=digits, na.print = "",
                  print.gap = 2)

    if(length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
        ## df component added in 1.8.0
        ## partial matching problem here.
        df <- if ("df" %in% names(x)) x[["df"]] else NULL
        if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
            cat("\nCoefficients: (", nsingular,
                " not defined because of singularities)\n", sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if(!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4L,
                            dimnames=list(cn, colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                     na.print="NA", ...)
    }
    ##
    cat("\n(Dispersion parameter for ", x$family$family,
	" family taken to be ", format(x$dispersion), ")\n\n","\nThe following statistics are valid only if this model is the FULL model in a model selection procedure\n",
	apply(cbind(paste(format(c("Null","Residual"), justify="right"),
                          "deviance:"),
		    format(unlist(x[c("null.deviance","deviance")]),
			   digits= max(5, digits+1)), " on",
		    format(unlist(x[c("df.null","df.residual")])),
		    " degrees of freedom\n"),
	      1L, paste, collapse=" "), sep="")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    cat("AIC: ", format(x$aic, digits= max(4, digits+1)),"\n\n",
	"Number of Fisher Scoring iterations: ", x$iter,
	"\n", sep="")

    correl <- x$correlation
    if(!is.null(correl)) {
# looks most sensible not to give NAs for undefined coefficients
#         if(!is.null(aliased) && any(aliased)) {
#             nc <- length(aliased)
#             correl <- matrix(NA, nc, nc, dimnames = list(cn, cn))
#             correl[!aliased, !aliased] <- x$correl
#         }
	p <- NCOL(correl)
	if(p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
		print(symnum(correl, abbr.colnames = NULL))
	    } else {
		correl <- format(round(correl, 2), nsmall = 2, digits = digits)
		correl[!lower.tri(correl)] <- ""
		print(correl[-1, -p, drop=FALSE], quote = FALSE)
	    }
	}
    }
    cat("\n")
    invisible(x)
}


