#############################################################
#
#   clm.* and related functions
#   Author: Claudio Agostinelli and Alessandro Gagliardi
#   E-mail: claudio@unive.it
#   Date: October, 30, 2013
#   Version: 0.1
#
#   Copyright (C) 2013 Claudio Agostinelli and Alessandro Gagliardi
#
#############################################################

clm.control <- function(epsilon = 1e-8, maxit = 50, nstart=2, trace = FALSE) {
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(nstart) || nstart != round(nstart))
    stop("'nstart' must be an integer")
  list(epsilon = epsilon, maxit = maxit, nstart = nstart, trace = trace)
}

clm <- function(formula, mu.est=TRUE, family = vonMises, data, weights,
		subset, na.action, center=FALSE, start = NULL,  
		mustart=NULL, kappastart=NULL, offset,
		control = list(...), control.circular = list(),
                model = TRUE, method = "clm.fit",
                x = FALSE, y = TRUE, contrasts = NULL, ...) {
  call <- match.call()
  ## family
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  ## extract x, y, etc from the model formula and frame
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
                  "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  if(identical(method, "model.frame"))
    return(mf)

  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  ## for back-compatibility in return result
  if (identical(method, "clm.fit"))
    control <- do.call("clm.control", control)

  mt0 <- mt <- attr(mf, "terms") # allow model.frame to have updated it
  attr(mt0,"intercept") <- 0
  if (is.empty.model(mt0) & mu.est)
    mt <- mt0
  Y <- model.response(mf, "any") # e.g. factors are allowed
  ## avoid problems with 1D arrays, but keep names
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  ## null model support
  X <- if (is.empty.model(mt))
	 matrix(,NROW(Y), 0L)
       else
         model.matrix(mt, mf, contrasts)
  
  posint <- which(colnames(X) == "(Intercept)")
  ## Centering the X?
  if (is.logical(center) && center)
    center <- apply(X[,-posint,drop=FALSE], 2, median)
  if (is.function(center))
    center <- apply(X[,-posint,drop=FALSE], 2, center)
  X[,-posint] <- scale(X[,-posint,drop=FALSE], center=center, scale=FALSE)

  ## avoid any problems with 1D or nx1 arrays by as.vector.
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  ## check weights and offset
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")

  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
  }
  ## these allow starting values to be expressed in terms of other vars.
  mustart <- model.extract(mf, "mustart")
  kappastart <- model.extract(mf, "kappastart")

  fit <- eval(call(if(is.function(method)) "method" else method,
    x = X, y = Y, weights = weights, start = start,
    mustart = mustart, kappastart = kappastart, offset = offset,
    family = family, control = control, mu.est = mu.est))

  if (length(offset) && mu.est) {
    fit2 <- eval(call(if(is.function(method)) "method" else method,
      x = matrix(,NROW(Y), 0L), y = Y, weights = weights, offset = offset,
      family = family, control = control, mu.est = mu.est))
    if (!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }
  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if(x)
    fit$x <- X
  if(!y)
    fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula, terms = mt,
    data = data, offset = offset, control = control, method = method,
    contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)))
  class(fit) <- c(fit$class, c("clm","glm", "lm"))
  fit
}

clm.fit <- function(x, y, mu.est=TRUE, weights = rep(1, nobs),
  start = NULL, mustart=NULL, kappastart=NULL, offset = rep(0, nobs),
  family = vonMises(), control = list(), control.circular = list()) {

  if (!is.logical(mu.est))
    stop("'mu.est' must be logical")
  
  control <- do.call("clm.control", control)

### Dependent variables
  if (is.circular(y))
    datacircularp <- circularp(y)
  else
    datacircularp <- list(type = "angles", units = "radians", template = "none",
                          modulo = "asis", zero = 0, rotation = "counter")
  dc <- control.circular
  if (is.null(dc$type)) 
    dc$type <- datacircularp$type
  if (is.null(dc$units)) 
    dc$units <- datacircularp$units
  if (is.null(dc$template)) 
    dc$template <- datacircularp$template
  if (is.null(dc$modulo)) 
    dc$modulo <- datacircularp$modulo
  if (is.null(dc$zero)) 
    dc$zero <- datacircularp$zero
  if (is.null(dc$rotation)) 
    dc$rotation <- datacircularp$rotation
  
  #4 here
  y <- conversion.circular(y, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
  attr(y, "circularp") <- attr(y, "class") <- NULL
  ynames <- names(y)
  nobs <- NROW(y)
## Explanatory variables  
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  conv <- FALSE
  nvars <- ncol(x)
  EMPTY <- nvars == 0
## define weights and offset if needed
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
## get family functions:
  variance <- family$variance
  linkinv  <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv) )
    stop("'family' argument seems not to be a valid family object", call. = FALSE)
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null) if(is.null(x)) if.null else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu  <- unless.null(family$validmu,  function(mu) TRUE)
  if(!is.null(mustart) & mu.est) {
    mustart <- conversion.circular(mustart, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
    attr(mustart, "circularp") <- attr(mustart, "class") <- NULL 
  }

  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta))
      stop("invalid linear predictor values in empty model", call. = FALSE)
    mulinear <- linkinv(eta)
    ## calculate initial deviance and coefficient
    if (!validmu(mulinear))
      stop("invalid fitted means in empty model", call. = FALSE)
    S <- sum(weights * sin(y-mulinear)) / sum(weights)
    C <- sum(weights * cos(y-mulinear)) / sum(weights)
    R <- sqrt((S^2+C^2))
    khat <- A1inv(R)
    if (mu.est)
      muhat <- atan2(S/R,C/R)
    else
      muhat <- 0
    dev <- sum(family$dev.resids(y, muhat, mulinear, khat, weights))
    
    S0 <- sum(weights * sin(y - linkinv(offset)))/sum(weights)
    C0 <- sum(weights * cos(y - linkinv(offset)))/sum(weights)
    R0 <- sqrt((S0^2+C0^2))
    khat0 <- A1inv(R0)
    muhat0 <- atan2(S0/R0,C0/R0)
      
    null.dev <- sum(family$dev.resids(y, muhat0, offset, khat0, weights))
    df.residual <- nobs - as.numeric(mu.est)
    result <- list(coefficients = numeric(), mu = muhat, kappa = khat,
      fitted.values = rep(muhat - mulinear, nobs),
      residuals = (y - muhat - mulinear), mu.est=mu.est,
      iter = 0L, family = family, aic = family$aic(dev=dev,rank=df.residual),
      deviance = dev, null.deviance = null.dev, y = y, R = R, rank = 0,
      df.residual = df.residual, df.null = nobs - as.numeric(mu.est),
      converged = TRUE, weights = weights, prior.weights = weights)
  } else {
    result <- ClmLocationCircularRad(x=x,y=y,mu.est=mu.est,weights=weights,offset=offset,start=start,mustart=mustart, kappastart=kappastart,family=family,epsilon=control$epsilon,maxit=control$maxit,trace=control$trace,nstart=control$nstart)
  }
    
  class(result) <- "clm"
  result$fitted.values <- conversion.circular(circular(result$fitted.values), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  result$mu <- minusPiPlusPi(conversion.circular(circular(result$mu), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation))
  result$residuals <- conversion.circular(circular(result$residuals), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  result$y <- conversion.circular(circular(result$y), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  attr(result$mu,"class") <- attr(result$residuals,"class") <- "circular"  
  return(result)
}

ClmLocationCircularRad <- function(x, y, mu.est=TRUE, weights, offset, start, mustart, kappastart, family, epsilon, maxit, trace, nstart) {
      
  logLikelihood <- function(khat,muhat,mulinear,weights,y) {
    if (khat < 100000)
      llik <- -sum(weights)*(log(2*pi)+log(besselI(khat,nu=0,expon.scaled=TRUE))+khat) + sum(weights * khat * cos(y-muhat-mulinear))
    else
      llik <- ifelse(all((y-muhat-mulinear)==0), sqrt(.Machine$double.xmax), -sqrt(.Machine$double.xmax))
    return(llik)
  }

  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  if (is.null(start) && nstart == 0L)
    stop("'nstart' must be positive if 'start' is null")

  ### initial value for beta
  beta <- ClmLocationCircularStartRad(x=x,y=y,mu.est=mu.est,weights=weights,offset=offset,start=start,mustart=mustart,kappastart=kappastart,family=family,nstart=nstart)
  
  lastbeta <- beta + 1 + epsilon
  nobs <- length(y)
  conv <- FALSE
  lll <- 0
  for (iter in 1L:maxit) {
#### Get muhat and khat
    eta <- drop(x%*%beta) + offset
    mulinear <- linkinv(eta)
    S <- sum(weights * sin(y - mulinear))/sum(weights)
    C <- sum(weights * cos(y - mulinear))/sum(weights)
    R <- sqrt((S^2+C^2))
    if (mu.est)
      muhat <- atan2(S/R,C/R)
    else
      muhat <- 0

    khat <- A1inv(R)
#### Get g values
    g <- mu.eta(eta)
#### Get u values
    uvector <- sin(y-muhat-mulinear)
#### Update Step
    lastbeta <- beta
    ustar <- uvector/(R*g)  ### R=A1(khat)
    fit <- lm.wfit(y=ustar,x=x,w=weights*(g^2))
    beta <- fit$coefficients + beta
    xgu <- t(x)%*%diag(weights)%*%diag(g)%*%uvector
    if (khat>=1e5) {
      khat<- 1e5-1
      conv <- TRUE
      break
    }
    if (any(is.na(beta))) {
      conv <- FALSE
      break
    }
    if (max(abs(xgu)) < 1e-7 & max(abs(beta-lastbeta)) > epsilon)
      beta <- (beta + lastbeta)/2
    if (trace) {
      cat("At iteration ",iter," :","\n")
      cat("log likelihood ",logLikelihood(khat=khat,muhat=muhat,mulinear=mulinear,weights=weights,y=y),"\n")
      cat("coefficients =", beta, "mu =", muhat, " kappa =", khat,"\n")
      cat("Value of equation=",t(x)%*%diag(weights)%*%diag(g)%*%uvector,"\n\n")
    }

    if (max(abs(beta-lastbeta)) < epsilon) {
      conv <- TRUE
      break
    } else
      conv <- FALSE
  }
  df.residual <- nobs - fit$rank - as.numeric(mu.est)
  df.null <- df.residual + fit$rank
  aic <- NA
  if (!conv & maxit > 1)
    warning("clm.fit: algorithm did not converge", call. = FALSE)
  if (khat==1e5-1)
    warning("clm.fit: Max 'kappa' value was reached", call. = FALSE)

  fitted.values <- muhat + linkinv(drop(x%*%beta)+offset)
  residuals <- y - fitted.values

  deviance <- sum(family$dev.resids(y=y, mu=muhat, mulinear=mulinear, kappa=khat, wt=weights))
  aic <- family$aic(dev=deviance,rank=df.residual)
  if (is.infinite(aic))
    aic <- -Inf

  S0 <- sum(weights * sin(y))/sum(weights)
  C0 <- sum(weights * cos(y))/sum(weights)
  R0 <- sqrt((S0^2+C0^2))
  khat0 <- A1inv(R0)
  muhat0 <- atan2(S0/R0,C0/R0)
  
  null.deviance <- sum(family$dev.resids(y=y, mu=muhat0, mulinear=0, kappa=khat0, wt=weights))
  
  result <- list(coefficients=drop(beta), mu=muhat, kappa=khat, fitted.values=fitted.values, residuals=residuals, mu.est=mu.est, iter=iter, family=family, aic=aic, deviance=deviance, null.deviance=null.deviance, y=y, R=R, rank = fit$rank, df.residual = df.residual, df.null = df.null, converged = conv, weights = weights*g^2, prior.weights = weights, effects=fit$effects, qr = fit$qr, linear.predictors = eta, xgu=xgu)
  return(result)
}

ClmLocationCircularStartRad <- function(x, y, mu.est, weights, offset, start, mustart, kappastart, family, nstart) {
      
  logLikelihood <- function(khat,muhat,mulinear,weights,y) {
    if (khat < 100000)
      llik <- -sum(weights)*(log(2*pi)+log(besselI(khat,nu=0,expon.scaled=TRUE))+khat) + sum(weights * khat * cos(y-muhat-mulinear))
    else
      llik <- ifelse(all((y-muhat-mulinear)==0), sqrt(.Machine$double.xmax), -sqrt(.Machine$double.xmax))
    return(llik)
  }

  beta.init <- function(x, y, mu, family, weights, offset, cc=2, ...) {
    z <- MinusPiPlusPiRad(y-mu)
    zz <- family$linkfun(z)
    az <- abs(z)
    w <- rep(0, length(z))
    pos <- which(az <=cc/2)
    w[pos] <- 1
    pos <- which(az > cc/2 & az <=cc)
    w[pos] <- (az[pos]-cc)^2
    w <- w*weights
    coeff <- lm(zz~x-1, weights=w, offset=offset, ...)$coeff
    return(coeff)
  }

  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  if (is.null(beta) && nstart == 0L)
    stop("'nstart' must be positive if 'start' is null")

  res <- MlevonmisesRad(y)
  ## if (!is.finite(res[4L]))
  ##   res[4L] <- 1e5-1L
  if (nstart > 2) {
    stepstart <- pi/(nstart+1)
    sequence <- seq(from=stepstart, to=pi - stepstart, by=stepstart)
    muinit <- res[1L]  +  c(0,sequence, -sequence, pi)	
  } else
    muinit <- c(0,pi) + res[1L]
  
  muinit <- MinusPiPlusPiRad(muinit)
  opt2 <- list(value=Inf)
  betainit <- matrix(0,nrow=length(muinit), ncol=NCOL(x))
  kappainit <- rep(0, length.out=length(muinit))
  for (i in 1:length(muinit)) {
    betainit[i,] <- beta.init(x=x, y=y, mu=muinit[i], family=family, weights=weights, offset=offset, cc=2)
    if (any(is.na(betainit[i,])))
      betainit[i,] <- rep(0, NCOL(x))
    eta <- drop(x%*%betainit[i,]) + offset
    mulinear <- linkinv(eta)
    res <- MlevonmisesRad(x=y - mulinear, mu=muinit[i])
    kappainit[i] <- res[4L]
  }
  betainit <- cbind(betainit, muinit, kappainit)
  if (!is.null(start)) {
    eta <- drop(x%*%start) + offset
    mulinear <- linkinv(eta)
  } else {
    mulinear <- 0
    start <- rep(0, NCOL(x))
  }
    
  res <- MlevonmisesRad(x=y - mulinear, mu=mustart, kappa=kappastart)
  mustart <- res[1L]
  kappastart <- res[4L]
  betainit <- rbind(betainit, c(start,mustart,kappastart))
  for (i in 1L:nrow(betainit)) {
    opt1 <- tryCatch(optim(par=betainit[i,], fn=function(x,y,z,wt) {
      nc <- NCOL(z)
      mulinear <- linkinv(drop(z%*%(x[1L:nc])))
      muhat <- x[nc + 1]
      khat <- x[nc + 2]
      llik <- -logLikelihood(khat = khat, muhat = muhat, mulinear=mulinear, weights=wt, y=y)
      return(llik)
    }, y=y, z=x, wt=weights, lower=c(rep(-Inf,NCOL(x)), -pi, 0.1), upper=c(rep(Inf,NCOL(x)), pi, Inf), method="L-BFGS-B"), error=function(e)
{
        return(list(value=0, par=rep(0,NCOL(x))))
})
    if (opt1$value < opt2$value) {
      beta <- drop(opt1$par[1L:NCOL(x)])
    } else {
      beta <- drop(opt2$par[1L:NCOL(x)])
    }
    opt2 <- opt1
  }
  return(beta)
}

print.clm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {  
  cat("\nCall:  ",
    paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x)) | x$mu.est) {
    cat("Coefficients")
    if (is.character(co <- x$contrasts))
      cat("  [contrasts: ",
        apply(cbind(names(co),co), 1L, paste, collapse = "="), "]")
    cat(":\n")
    coeff <- x$coefficients
    if (x$mu.est) {
      class(x$mu) <- "numeric"
      attr(x$mu, "circularp") <- NULL
      name <- names(coeff)
      coeff <- c(x$mu, coeff)
      names(coeff) <- c("(Mu Est.)",name)
    }
    print.default(format(coeff, digits = digits), print.gap = 2, quote = FALSE)

  } else
    cat("No coefficients\n\n")

  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
    x$df.residual, "Residual\n")
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (",mess, ")\n", sep = "")
  cat("AIC:", format(signif(x$aic, digits)), "\n")
  invisible(x)
}

summary.clm <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...) {
  df.r <- object$df.residual
  dispersion <- 1/(object$kappa*object$R)
  ## calculate scaled and unscaled covariance matrix
  intercept <- attr(object$terms, "intercept")

  p <- object$rank
  aliased <- is.na(coef(object)) # used in print method
  if (object$mu.est) {
        aliased <- c(aliased, is.na(object$mu))
        anames <- names(aliased)
        names(aliased)[length(aliased)] <- "(Mu Est.)"
  }

  if (p > 0) {
    p1 <- 1L:p
    qr.clm <- function(x, ...) {
      if (is.null(x$qr))
        stop("clm object does not have a proper 'qr' component.\nRank zero or should not have used clm(.., qr=FALSE).")
      x$qr
    }
    Qr <- qr.clm(object)
    coef.p <- object$coefficients[Qr$pivot[p1]]
    covmat.unscaled <- chol2inv(Qr$qr[p1,p1,drop=FALSE])
    dimnames(covmat.unscaled) <- list(names(coef.p),names(coef.p))
    covmat <- dispersion*covmat.unscaled
    var.cf <- diag(covmat)
    ## calculate coef table

    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    if (object$mu.est) {
      s.mu.err <- sqrt(dispersion/(object$df.residual+2))
      if (attr(object$mu, "circularp")$units=="degrees")
	s.mu.err <- s.mu.err * 180/pi
      mu.tvalue <- object$mu/s.mu.err
    }
    dn <- c("Estimate", "Std. Error")
    rn <- names(coef.p)
    if (df.r > 0) {
      pvalue <- 2*pt(-abs(tvalue), df.r)
      coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
      if (object$mu.est) {
        mu.pvalue <- 2 * pt(-abs(mu.tvalue), df.r)				
        mu.table <- c(as.numeric(object$mu), s.mu.err, mu.tvalue, mu.pvalue)
        coef.table <- rbind(mu.table, coef.table)
        rn <- c("(Mu Est.)", rn)
      }		
      dimnames(coef.table) <- list(rn, c(dn, "t value","Pr(>|t|)"))
    } else { # df.r == 0
      if (object$mu.est) {
        coef.p <- c(object$mu, coef.p)
        rn <- c("(Mu Est.)", rn)
      }
      coef.table <- cbind(coef.p, NaN, NaN, NaN)
      dimnames(coef.table) <- list(rn, c(dn, "t value","Pr(>|t|)"))
    }
    df.f <- NCOL(Qr$qr)
  } else {
    coef.table <- matrix(,as.numeric(object$mu.est), 4L)	
    rn <- NULL
    covmat.unscaled <- covmat <- matrix(, 0L, 0L)
    if (object$mu.est) {
      s.mu.err <- sqrt(dispersion/(object$df.residual+2))
      if (attr(object$mu, "circularp")$units=="degrees")
        s.mu.err <- s.mu.err * 180/pi
      mu.tvalue <- object$mu/s.mu.err
      mu.pvalue <- 2 * pt(-abs(mu.tvalue), df.r)
      coef.table[1,] <- c(as.numeric(object$mu),s.mu.err,mu.tvalue,mu.pvalue)
      rn <- "(Mu Est.)"
    }
    dimnames(coef.table) <-
      list(rn, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    df.f <- 0
  }
  ## return answer

  ## these need not all exist, e.g. na.action.
  keep <- match(c("call","terms","family","deviance", "aic",
		    "contrasts", "df.residual","null.deviance","df.null",
                    "iter", "na.action"), names(object), 0L)
    
  ans <- c(object[keep],
    list(deviance.resid = residuals(object["residuals"], type = "deviance"),
    coefficients = coef.table,
    aliased = aliased,
    dispersion = dispersion,
    kappa=object$kappa,
    df = c(object$rank, df.r, df.f),
    cov.unscaled = covmat.unscaled,
    cov.scaled = covmat))

  if (correlation && p > 0) {
    dd <- sqrt(diag(covmat.unscaled))
    ans$correlation <-
      covmat.unscaled/outer(dd,dd)
    ans$symbolic.cor <- symbolic.cor
  }
  class(ans) <- "summary.clm"
  return(ans)
}

print.summary.clm <- function (x, digits = max(3L, getOption("digits") - 3L),
	      symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Residuals: \n")
  if (x$df.residual > 5 & sum(is.na(x$deviance.resid)) < 5) {
    summre <- as.numeric(minusPiPlusPi(quantile(x$deviance.resid,na.rm = TRUE)))
    summre <- setNames(summre, c("Min", "1Q", "Median", "3Q", "Max"))
  }
  else
  {
    summre <- as.numeric(quantile(x$deviance.resid,na.rm = TRUE))
  }
  xx <- zapsmall(summre, digits + 1L)
  print.default(xx, digits = digits, na.print = "", print.gap = 2L)

  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  } else {
## partial matching problem here.
    df <- if ("df" %in% names(x)) x[["df"]] else NULL
    if ((!is.null(df) && (nsingular <- df[3L] - df[1L])))
      cat("\nCoefficients: (", nsingular,
        " not defined because of singularities)\n", sep = "")
    else
      cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L,
        dimnames=list(cn, colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
      na.print = "NA", ...)
  }
##
  cat("\nConcentration parameter for ", x$family$family,
    " family is ", format(x$kappa), "\n\n",
    apply(cbind(paste(format(c("Null","Residual"), justify="right"),
    "deviance:"), format(unlist(x[c("null.deviance","deviance")]),
    digits = max(5L, digits + 1L)), " on",
    format(unlist(x[c("df.null","df.residual")])),
    " degrees of freedom\n"), 1L, paste, collapse = " "), sep = "")
  if(nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),"\n\n",
    "Number of Fisher Scoring iterations: ", x$iter,
    "\n", sep = "")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      } else {
        correl <- format(round(correl, 2L), nsmall = 2L, digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop=FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}

model.frame.clm <- function (formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    fcall$method <- "model.frame"
    fcall[[1L]] <- as.name("clm")
    fcall[names(nargs)] <- nargs
##  env <- environment(fcall$formula)  # always NULL
    env <- environment(formula$terms)
    if (is.null(env))
      env <- parent.frame()
    eval(fcall, env)
  } else formula$model
}

weights.clm <- function(object, type = c("prior", "working"), ...) {
  type <- match.arg(type)
  res <- if(type == "prior") object$prior.weights else object$weights
  if(is.null(object$na.action))
    res
  else
    naresid(object$na.action, res)
}

formula.clm <- function(x, ...) {
  form <- x$formula
  if ( !is.null(form) ) {
    form <- formula(x$terms) # has . expanded
    environment(form) <- environment(x$formula)
    form
  } else formula(x$terms)
}
