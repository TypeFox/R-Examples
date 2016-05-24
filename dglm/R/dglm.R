dglm <- function(formula = formula(data),
                 dformula = ~1,
                 family = stats::gaussian,
                 dlink = "log",
                 data = sys.parent(),
                 subset = NULL,
                 weights = NULL,
                 contrasts = NULL,
                 method = "ml",
                 mustart = NULL,
                 betastart = NULL,
                 etastart = NULL,
                 phistart = NULL,
                 control = dglm.control(...),
                 ykeep = TRUE,
                 xkeep = FALSE,
                 zkeep = FALSE,
                 ...) {
  #   Double generalized linear models
  #   Gordon Smyth, Walter and Eliza Hall Institute of Medical Research
  #   S-Plus version created 8 Dec 1997, last revised 22 Oct 1999.
  #
  #   Ported to R by Peter Dunn, 22 August 2005
  #
  #  Set up mean submodel: 
  #               y   response
  #               X   design matrix
  #               w   prior weights
  #          offset   offset in linear predictor
  #          family   response family
  #         mustart   starting values for mu (optional)
  #       betastart   starting values for coefficients (optional)
  #
  #   Save call for future reference
  call <- match.call()
  
  #   Get family for mean model
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  #   Evaluate the model frame
  mnames <- c("", "formula", "data", "weights", "subset")
  cnames <- names(call)
  
  cnames <- cnames[match(mnames,cnames,0)]
  mcall <- call[cnames]
  mcall[[1]] <- as.name("model.frame")
  mean.mframe <- eval(mcall, parent.frame())  
    ### Dunn:  NEEDS THE <<- ############
    ### Corty: I copied this line from glm() and it seems to work
    ###        See notes in git repo.  Solved a problem w R CMD CHECK
  # mf <- match.call(expand.dots = FALSE)  # commented this out with seemingly no problem.
  
  #   Now extract the glm components
  y <- stats::model.response(mean.mframe, "numeric")
  
  if (is.null(dim(y))) {
    N <- length(y) 
  } else {
    N <- dim(y)[1]
  }
  nobs <- N # Needed for some of the  eval  calls
  mterms <- attr(mean.mframe, "terms")
  X <- stats::model.matrix(mterms, mean.mframe, contrasts)
  
  weights <- stats::model.weights(mean.mframe)
  if (is.null(weights)) weights <- rep(1,N)
  if (!is.null(weights) && any(weights < 0)) {
    stop("negative weights not allowed")
  }
  
  offset <- stats::model.offset(mean.mframe)
  if ( is.null(offset) ) offset <- rep(0, N )
  if (!is.null(offset) && length(offset) != NROW(y)) {
    stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                  length(offset), NROW(y)), domain = NA)
  }
  
  #
  #  Set up dispersion submodel:
  #               Z   design matrix 
  #         doffset   offset in linear predictor
  #         dfamily   family for unit deviances
  #        phistart   starting values for phi (optional)
  #
  #   Setup dformula with y on left hand side
  #   (to ensure that mean and dispersion submodels have same length)
  mcall$formula <- formula
  mcall$formula[3] <- switch(match(length(dformula),c(0,2,3)),
                             1,dformula[2],dformula[3])
  
  #   Evaluate the model frame and extract components
  var.mframe <- eval(mcall, sys.parent())
  dterms <- attr(var.mframe, "terms")
  
  Z <- stats::model.matrix(dterms, var.mframe, contrasts)
  doffset <- stats::model.extract(var.mframe, offset)
  if ( is.null(doffset) ) doffset <- rep(0, N )

  #   Parse the dispersion link and evaluate as list
  name.dlink <- substitute(dlink)
  if(is.name(name.dlink)) {
    if(is.character(dlink)) { # link is character variable
      name.dlink <- dlink
    } else {                  # link is name without quotes
      dlink <- name.dlink <- as.character(name.dlink)
    }
  } else {
    if(is.call(name.dlink))  # power link
      name.dlink <- deparse(name.dlink)
  }
  if(!is.null(name.dlink))   name.dlink <- name.dlink

  #   Construct the dispersion variance function
  if ( family$family=="Tweedie") {
    tweedie.p <- call$family$var.power
  }
  
  Digamma <- family$family=="Gamma" || (family$family=="Tweedie" && tweedie.p==2) 
  # In other words, if the mean glm is a gamma (or, equivalently, Tweedie with p=2), use the digamma explicitly
  if (Digamma) {
    
    linkinv <- stats::make.link(name.dlink)$linkinv
    linkfun <- stats::make.link(name.dlink)$linkfun
    mu.eta <- stats::make.link(name.dlink)$mu.eta
    valid.eta <- stats::make.link(name.dlink)$valid.eta
    init <- expression({
      if (any(y <= 0)){
        print(y)
        print( any(y<=0) )
        stop("non-positive values not allowed for the DM gamma family")
      }
      n <- rep.int(1, nobs)
      mustart <- y
    })
    
    
    dfamily <- structure(list(family="Digamma", variance=varfun.digamma, 
                              dev.resids=function(y,mu,wt){wt * unitdeviance.digamma(y,mu)}, # gaussian()$dev.resids are resids^2
                              aic=function(y, n, mu, wt, dev) NA,  
                              link=name.dlink, linkfun=linkfun, linkinv=linkinv, 
                              mu.eta=mu.eta, initialize=init, 
                              validmu=function(mu) { all(mu>0) }, valideta=valid.eta) )
    
  } else { # Not digamma family
    eval(substitute(dfamily <- Gamma(link=lk), list(lk=name.dlink ) ))
  }
  #
  #   Remember if log-link for use with initial values
  #
  dlink <- as.character(dfamily$link)
  logdlink <- dlink=="log"
  #
  #  Match method (ml or reml)
  #
  if(!is.null(call$method)) {
    name.method <- substitute(method)
    if(!is.character(name.method))
      name.method <- deparse(name.method)
    list.methods <- c("ml","reml","ML","REML","Ml","Reml")
    i.method <- pmatch(method,list.methods,nomatch=0)
    if(!i.method) stop("Method must be ml or reml")
    method <- switch(i.method,"ml","reml","ml","reml","ml","reml")
  }
  reml <- method=="reml"
  #
  #  Starting values.  If explicit starting values are not supplied,
  #  regress link(y) on X and dlink(d) on Z by ordinary least squares.
  
  if ( is.null(mustart) ) { 
    etastart <- NULL
    eval(family$initialize)  
    mu <- mustart
    mustart <- NULL
  }
  
  if(!is.null(betastart)) {
    eta <- X %*% betastart
    mu <- family$linkinv(eta+offset)
  } else {
    if(!is.null(mustart)) {
      mu <- mustart
      eta <- family$linkfun(mu)-offset
    } else {
      
      eta <- stats::lm.fit(X,family$linkfun(mu)-offset,singular.ok=TRUE)$fitted.values
      # Recall:  fitted values are on the linear predictor scale
      mu <- family$linkinv(eta+offset)
    }
  }
  
  #   Careful:   Called dev.resid, but appear to be the deviance residuals squared
  d <- family$dev.resids(y, mu, weights)
  
  if (!is.null(phistart)) {
    phi <- phistart
    deta <- dfamily$linkfun(phi) - doffset
  } else {
    deta <- stats::lm.fit(Z,dfamily$linkfun(d + (d == 0)/6) - doffset,singular.ok = TRUE)$fitted.values
    if (logdlink) deta <- deta + 1.27036
    phi <- dfamily$linkinv(deta + offset)
  }
  if (any(phi <= 0)) {
    cat("Some values for  phi  are non-positive, suggesting an inappropriate model",
        "Try a different link function.\n")
  }
  
  zm <- as.vector( eta + (y - mu) / family$mu.eta(eta) )
  wm <- as.vector( eval(family$variance(mu))*weights/phi )
  # as.vector()  added 30 October 2012
  mfit <- stats::lm.wfit(X, zm, wm, method = "qr", singular.ok = TRUE) # ,qr=reml)
  eta <- mfit$fitted.values
  
  
  mu <- family$linkinv(eta + offset)
  if ( family$family == "Tweedie") {
    cat("p:",tweedie.p,"\n")
    if ( (tweedie.p > 0) & (any(mu < 0)) ) {
      cat("Some values for  mu  are negative, suggesting an inappropriate model.",
          "Try a different link function.\n")
    }
  }
  d <- family$dev.resids(y, mu, weights)
  
  #  Initial (minus twice log) likelihood or adjusted profile likelihood
  const <- dglm.constant(y,family,weights)
  
  if (Digamma) {
    h <- 2*(lgamma(weights/phi) + (1 + log(phi/weights))*weights/phi)
  } else {
    h <- log(phi/weights)
  }
  m2loglik <- const + sum(h + d/phi)
  
  if (reml)
    m2loglik <- m2loglik + 2*log(abs(prod(diag(mfit$R))))
  
  m2loglikold <- m2loglik + 1
  #  Estimate model by alternate iterations
  
  epsilon <- control$epsilon
  maxit <- control$maxit
  trace <- control$trace
  iter <- 0
  
  while (abs(m2loglikold-m2loglik)/(abs(m2loglikold) + 1) > epsilon && iter < maxit )  {
    
    ################################
    #      dispersion submodel
    hdot <- 1/dfamily$mu.eta( deta )
    
    if (Digamma) {
      delta <- 2*weights*(log(weights/phi) - digamma(weights/phi))
      u <- 2*(weights ^ 2)*(trigamma(weights/phi) - phi/weights)
      fdot <- (phi ^ 2) / u * hdot  # p 50
    } else {# In normal and iG cases, eg, the dispersion sub-model is gamma
      delta <- phi
      u <- (phi ^ 2) # variance function for disp. model; u(delta)=delta^2 is gamma
      fdot <- hdot
    }
    wd <- 1 / ((fdot ^ 2) * u)   # ie Smyth, p 50.  We don't include the factor of 2,
    # as that is the disp. parameter, which enters via
    # the  dispersion=2  argument for the summary.
    
    if(reml) {
      h <- stats::hat(mfit$qr)
      delta <- delta - phi*h
      wd <- wd - 2*(h/hdot^2/phi^2) + h^2
    }
    
    if(any(wd<0)) {
      cat(" Some weights are negative; temporarily fixing.  This may be a sign of an inappropriate model.\n")
      wd[wd<0] <- 0
    }
    if(any(is.infinite(wd))) {
      cat(" Some weights are negative; temporarily fixing.  This may be a sign of an inappropriate model.\n")
      wd[is.infinite(wd)] <- 100
    }
    zd <- deta + (d - delta) * fdot
    
    # Now fit dispersion submodel, with response zd, explanatory vars Z, weights wd
    dfit <- stats::lm.wfit(Z, zd, wd, method="qr", singular.ok=TRUE)
    deta <- dfit$fitted.values
    
    phi <- dfamily$linkinv(deta+doffset)
    if (any(is.infinite(phi))) {
      cat("*** Some values for  phi  are infinite, suggesting an inappropriate model",
          "Try a different link function.  Making an attempt to continue...\n")
      phi[is.infinite(phi)] <- 10
    }
    
    ################################
    #      mean submodel
    zm <- eta + (y - mu) / family$mu.eta(eta)
    fam.wt <- expression( weights * family$variance(mu) ) 
    wm <- eval( fam.wt )/phi
    mfit <- stats::lm.wfit(X, zm, wm, method = "qr", singular.ok = TRUE)
    eta <- mfit$fitted.values
    mu <- family$linkinv(eta+offset)
    if (family$family=="Tweedie"){
      if ( (tweedie.p > 0) & (any(mu<0)) ) {
        cat("*** Some values for  mu  are negative, suggesting an inappropriate model.",
            "Try a different link function.  Making an attempt to continue...\n")
        mu[mu<=0] <- 1
      }
    }
    d <- family$dev.resids(y, mu, weights)
    
    
    #      overall likelihood
    m2loglikold <- m2loglik
    if(Digamma) {
      h <- 2*(lgamma(weights/phi)+(1+log(phi/weights))*weights/phi)
    } else {
      h <- log(phi/weights)
    }
    m2loglik <- const + sum(h + d/phi)
    if(reml) {
      m2loglik <- m2loglik + 2*log(abs(prod(diag(mfit$R))))
    }
    iter <- iter+1
    if(trace)
      cat("DGLM iteration ", iter, ": -2*log-likelihood = ",
          format(round(m2loglik, 4)), " \n", sep = "")
  } ### END while LOOP
  
  #
  #  Output for mean model:
  #  Exactly as for glm.object.  As for lm.object except that
  #  linear.predictors and prior.weights are new components, and fitted.values
  #  has a new definition.
  #
  mfit$call <- call
  mfit$formula <- formula #environment(call$formula)
  mfit$terms <- mterms
  
  mfit$model <- mean.mframe
  
  mfit$family <- family
  
  mfit$linear.predictors <- mfit$fitted.values+offset
  mfit$fitted.values <- mu
  mfit$prior.weights <- weights
  mfit$contrasts <- attr(X, "contrasts")
  intercept <- attr(mterms, "intercept")
  mfit$df.null <- N - sum(weights == 0) - as.integer(intercept)
  mfit$deviance <- sum(d/phi)
  mfit$aic <- NA
  mfit$null.deviance <- stats::glm.fit(x = X, y = y, weights = weights/phi, offset = offset, family = family)
  if (length(mfit$null.deviance) > 1) mfit$null.deviance <- mfit$null.deviance$null.deviance
  if (ykeep) mfit$y <- y
  if (xkeep) mfit$x <- X
  class(mfit) <- c("glm","lm")
  
  #
  #  Output for dispersion model:
  #  As for glm.object except that prior.weights are not relevant.  Is
  #  nested in one output component.
  #
  
  call$formula <- dformula
  dfit$terms <- dterms
  dfit$model <- var.mframe
  dfit$family <- dfamily
  dfit$prior.weights <- rep(1, N)
  dfit$linear.predictors <- dfit$fitted.values + doffset
  dfit$fitted.values <- phi
  dfit$aic <- NA
  call$dformula <- NULL
  call$family <- call(dfamily$family,link = name.dlink)
  dfit$call <- call
  dfit$residuals <- dfamily$dev.resid(d, phi, wt = rep(1/2,N) )
  dfit$deviance <- sum( dfit$residuals  )
  dfit$null.deviance <- stats::glm.fit(x = Z, y = d, weights = rep(1/2, N), offset = doffset, family = dfamily)
  if (length(dfit$null.deviance) > 1) dfit$null.deviance <- dfit$null.deviance$null.deviance
  if (ykeep) dfit$y <- d
  if (zkeep) dfit$z <- Z
  dfit$formula <- as.vector(attr(dterms, "formula"))
  dfit$iter <- iter
  class(dfit) <- c("glm","lm")
  out <- c(mfit, list(dispersion.fit = dfit, iter = iter, method = method, m2loglik = m2loglik))
  class(out) <- c("dglm","glm","lm")
  out
}
