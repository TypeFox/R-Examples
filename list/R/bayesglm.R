n <- NULL

bayesglm.internal <- function (formula, family = gaussian, data, weights, subset,
                      na.action, start = NULL, etastart, mustart, offset, control = glm.control(...),
                      method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL,
                      drop.unused.levels = TRUE,
                      scaled = TRUE, keep.order = TRUE,
                      drop.baseline = TRUE, n.iter = 100, 
                      print.unnormalized.log.posterior = FALSE, Warning=TRUE,...)
{
  call <- match.call()
  if (is.character(family)){
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)){
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data)){
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- drop.unused.levels
  mf$na.action <- NULL
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  #switch(method, model.frame = return(mf), glm.fit = 1, stop("invalid 'method': ", method))
  if (identical(method, "model.frame")) 
    return(mf)
  if (!is.character(method) && !is.function(method)) 
    stop("invalid 'method' argument")
  if (identical(method, "glm.fit")) 
    control <- do.call("glm.control", control)
  
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)){
      names(Y) <- nm
    }
  }
  # 2012.10.3  I input data instead of mf here.  Don't know if this is right
  X <- if (!is.empty.model(mt)) {
    model.matrixBayes.internal(object=mt, data=data, contrasts.arg=contrasts, keep.order = keep.order, drop.baseline=drop.baseline)
    #model.matrix.default(mt, mf, contrasts)
  }
  else{
    matrix(, NROW(Y), 0)
  }
  #  nobs <- NROW(X)
  #  X <- rbind(X, diag(NCOL(X)))
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) 
    stop("negative weights not allowed")
  
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                    length(offset), NROW(Y)), domain = NA)
  }
  mustart <- model.extract(mf, "mustart")
  
  etastart <- model.extract(mf, "etastart")
  
  fit <- bayesglm.fit.internal(x = X, y = Y, weights = weights, start = start,
                      family = family, control = glm.control(maxit = n.iter),
                      scaled = scaled)
  
  fit$model <- mf
  
  fit$na.action <- attr(mf, "na.action")
  if (x){
    fit$x <- X
  }
  if (!y){
    fit$y <- NULL
  }
  fit <- c(fit, list(call = call, formula = formula, terms = mt,
                     data = data, offset = offset, control = control, method = method,
                     contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)),
           keep.order = keep.order, drop.baseline = drop.baseline
  )
  class(fit) <- c("bayesglm", "glm", "lm")
  return(fit)
}

bayesglm.fit.internal <- function (x, y, weights = rep(1, nobs), start = NULL, 
        etastart = NULL, mustart = NULL, offset = rep(0, nobs), 
        family = gaussian(),
        control = glm.control(), intercept = TRUE,
        prior.mean = 0,
        prior.scale = NULL,
        prior.df = 1,
        prior.mean.for.intercept = 0,
        prior.scale.for.intercept = NULL,
        prior.df.for.intercept = 1,
        min.prior.scale=1e-12,
        scaled = TRUE, print.unnormalized.log.posterior=FALSE, Warning=TRUE)
{
    ######### INITIALIZE #######
    nobs <- NROW(y)
    nvars <- NCOL(x)
    conv <- FALSE
    EMPTY <- nvars == 0
    
    output <- .bayesglm.fit.initialize.priors (family, prior.scale, prior.scale.for.intercept, nvars, 
            prior.mean, prior.mean.for.intercept, intercept, prior.df, prior.df.for.intercept)
    prior.scale <- output$prior.scale 
    prior.scale.for.intercept <- output$prior.scale.for.intercept
    prior.mean <- output$prior.mean
    prior.scale <- output$prior.scale
    prior.df <- output$prior.df
    
    prior.scale <- .bayesglm.fit.initialize.priorscale (scaled, family, prior.scale, y, nvars, x, min.prior.scale)
    
    output <- .bayesglm.fit.initialize.x (x, nvars, nobs, intercept, scaled)
    x <- output$x
    xnames <- output$xnames
    x.nobs <- output$x.nobs
    
    
    output <- .bayesglm.fit.initialize.other (nobs, y, weights, offset)
    ynames <- output$ynames
    weights <- output$weights
    offset <- output$offset
    
    family <- .bayesglm.fit.initialize.family (family)
    ## TODO: DL -- I would put this inside initialize.family, but this does some magic setting of variables.
    ##             What I can do is push an environment into the function and wrap it in.
    if (is.null(mustart)){
        eval(family$initialize)
    }
    else {
        mustart.keep <- mustart
        eval(family$initialize)
        mustart <- mustart.keep
    }
    
    
    
    ########################
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!family$valideta(eta)){
            stop("invalid linear predictor values in empty model")
        }
        mu <- family$linkinv(eta)
        if (!family$validmu(mu)){
            stop("invalid fitted means in empty model")
        }
        devold <- sum(family$dev.resids(y, mu, weights))
        w <- ((weights * family$mu.eta(eta)^2)/family$variance(mu))^0.5
        residuals <- (y - mu)/family$mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
    }
    else {
        ## 3.1 ##
        output <- ..bayesglm.fit.loop.setup (etastart, start, nvars, xnames, offset, x, nobs, family, 
                mustart, eta, weights, conv, prior.scale, y, x.nobs) 
        coefold <- output$coefold
        eta <- output$eta
        mu <- output$mu
        devold <- output$devold
        boundary <- output$boundary
        prior.sd <- output$prior.sd
        dispersion <- output$dispersion
        dispersionold <- output$dispersionold
        ##########
        
        ## 3.2 ## check the link between 3.1 and 3.2
        output.new <- .bayesglm.fit.loop.main.ideal (control, x, y, nvars, nobs, weights, offset,
                intercept, scaled, 
                start = coefold, etastart = eta, mustart = mu,
                family,
                prior.mean, prior.mean.for.intercept, prior.scale, prior.df, prior.df.for.intercept,
                print.unnormalized.log.posterior,
                Warning)
        fit <- output.new$fit
        good <- output.new$good
        z <- output.new$z
        w <- output.new$w
        ngoodobs <- output.new$ngoodobs
        prior.scale <- output.new$prior.scale
        prior.sd <- output.new$prior.sd
        eta <- output.new$eta
        mu <- output.new$mu
        dev <- output.new$dev
        dispersion <- output.new$dispersion
        start <- output.new$start
        coef <- output.new$coef
        conv <- output.new$conv
        iter <- output.new$iter
        boundary <- output.new$boundary
        Rmat <- output.new$Rmat
        residuals <- output.new$residuals
    }
    ## 4 ##
    output <- .bayesglm.fit.cleanup (ynames, residuals, mu, eta, nobs, weights, w, good, weights, family$linkinv, family$dev.resids, y, intercept, fit, offset, EMPTY, n, dev, family$aic)
    residuals <- output$residuals
    mu <- output$mu
    eta <- output$eta
    wt <- output$wt
    weights <- output$weights
    y <- output$y
    wtdmu <- output$wtdmu
    nulldev <- output$nulldev
    n.ok <- output$n.ok
    nulldf <- output$nulldf
    rank <- output$rank
    resdf <- output$resdf
    aic.model <- output$aic.model
    #######
    list(coefficients = coef, 
            residuals = residuals, 
            fitted.values = mu,
            effects = if (!EMPTY) fit$effects, 
            R = if (!EMPTY) Rmat,
            rank = rank, 
            qr = if (!EMPTY) structure(getQr(fit)[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"), 
            family = family,
            linear.predictors = eta, 
            deviance = dev, 
            aic = aic.model,
            null.deviance = nulldev, 
            iter = iter, 
            weights = wt, 
            prior.weights = weights,
            df.residual = resdf, 
            df.null = nulldf, 
            y = y, 
            converged = conv,
            boundary = boundary, 
            prior.mean = prior.mean, 
            prior.scale = prior.scale,
            prior.df = prior.df, 
            prior.sd = prior.sd, 
            dispersion = dispersion)
}




.bayesglm.fit.initialize.priors <- function (family, prior.scale, prior.scale.for.intercept, nvars, 
        prior.mean, prior.mean.for.intercept, intercept, prior.df, prior.df.for.intercept) {
    ##### 12.13 ####
    if (is.null(prior.scale)){
        prior.scale <- 2.5
        if (family$link == "probit"){
            prior.scale <- prior.scale*1.6
        }
    }
    
    if (is.null(prior.scale.for.intercept)){
        prior.scale.for.intercept <- 10
        if (family$link == "probit"){
            prior.scale.for.intercept <- prior.scale.for.intercept*1.6
        }
    }
    ################
    
    if (length(prior.mean) == 1) {
        prior.mean <- rep(prior.mean, nvars)
    }
    else {
        if (intercept) {
            prior.mean <- c(prior.mean.for.intercept, prior.mean)
        }
    }
    
    if (length(prior.scale)==1){
        prior.scale <- rep(prior.scale, nvars)
    }
    else {
        if (intercept) {
            prior.scale <- c(prior.scale.for.intercept, prior.scale)
        }
    }
    
    if (length(prior.df) == 1) {
        prior.df <- rep(prior.df, nvars)
    }
    else {
        if (intercept) {
            prior.df <- c(prior.df.for.intercept, prior.df)
        }
    }
    
    list (prior.scale=prior.scale, 
            prior.scale.for.intercept=prior.scale.for.intercept,
            prior.mean=prior.mean,
            prior.scale=prior.scale,
            prior.df=prior.df)
}

.bayesglm.fit.initialize.priorscale <- function (scaled, family, prior.scale, y, nvars, x, min.prior.scale) {
    if (scaled) {
        if (family$family == "gaussian"){
            prior.scale <- prior.scale * 2 * sd(y)
        }
        
        for (j in 1:nvars) {
            x.obs <- x[, j]
            x.obs <- x.obs[!is.na(x.obs)]
            num.categories <- length(unique(x.obs))
            x.scale <- 1
            if (num.categories == 2) {
                x.scale <- max(x.obs) - min(x.obs)
            }
            else if (num.categories > 2) {
                x.scale <- 2 * sd(x.obs)
            }
            prior.scale[j] <- prior.scale[j]/x.scale
            
            if (prior.scale[j] < min.prior.scale){
                prior.scale[j] <- min.prior.scale
                warning ("prior scale for variable ", j,
                        " set to min.prior.scale = ", min.prior.scale,"\n")
            }
        }
    }
    return (prior.scale)
}

.bayesglm.fit.initialize.x <- function (x, nvars, nobs, intercept, scaled) {
    x <- as.matrix (rbind(x, diag(nvars)))
    xnames <- dimnames(x)[[2]]
    x.nobs <- x[1:nobs, ,drop=FALSE]
    # TODO: **** I moved it here because I think it's right, and changed it to x.nobs
    if (intercept & scaled) {
        x[nobs+1,] <- colMeans(x.nobs)
    }
    
    list (x=x, xnames=xnames, x.nobs=x.nobs)
}

.bayesglm.fit.initialize.family <- function (family, mustart, enviroment) {
    if (!is.function(family$variance) || !is.function(family$linkinv)){
        stop("'family' argument seems not to be a valid family object")
    }
    if (is.null(family$valideta)){
        family$valideta <- function(eta) TRUE
    }
    
    if (is.null(family$validmu)){
        family$validmu <- function(mu) TRUE
    }
    
    return (family)
}

.bayesglm.fit.initialize.other <- function (nobs, y, weights, offset) {
    if (is.matrix(y)){
        ynames <- rownames(y)
    }
    else{
        ynames <- names(y)
    }
    
    if (is.null(weights)){
        weights <- rep.int(1, nobs)
    }
    
    if (is.null(offset)){
        offset <- rep.int(0, nobs)
    }
    
    list (ynames=ynames,
            weights=weights,
            offset=offset)
}

..bayesglm.fit.loop.setup <- function (etastart, start, nvars, xnames, offset, x, nobs, family, 
        mustart, eta, weights, conv, prior.scale, y, x.nobs) {
    coefold <- NULL
    if (!is.null(etastart)) {
        eta <- etastart
    }
    else if (!is.null(start)) {
        if (length(start) != nvars)
            stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                            nvars, paste(deparse(xnames), collapse = ", ")), domain = NA)
        else {
            eta <- offset + as.vector(ifelse((NCOL(x) == 1), x.nobs[,1]*start, x.nobs %*% start))
            coefold <- start
        }
    }
    else {
        eta <- family$linkfun(mustart)
    }
    
    mu <- family$linkinv(eta)
    if (!isTRUE(family$validmu(mu) && family$valideta(eta)))
        stop("cannot find valid starting values: please specify some")
    devold <- sum(family$dev.resids(y, mu, weights))
    boundary <- conv <- FALSE
    prior.sd <- prior.scale
    #========Andy 2008.7.8=============
    dispersion <- ifelse((family$family %in% c("poisson", "binomial")),  1, var(y)/10000)
    #==================================
    dispersionold <- dispersion
    
    list (coefold=coefold,
            eta=eta,
            mu=mu,
            devold=devold,
            boundary=boundary,
            prior.sd=prior.sd,
            dispersion=dispersion,
            dispersionold=dispersionold)
}


#.bayesglm.fit.loop.initialMemoryAllocation <- function (epsilon, nvars, ngoodobs) {
#    list (ny=as.integer(1),
#            tol = min(1e-07, epsilon/1000),
#            coefficients = double(nvars),
#            rank = integer(1),
#            pivot = 1:nvars,
#            qraux = double (nvars),
#            work = double (2 * nvars),
#            residuals = double (ngoodobs + nvars),
#            effects=double (ngoodobs + nvars))
#}

.bayesglm.fit.loop.initializeState <- function (start, etastart, mustart, offset, x.nobs, var.y, nvars, family, weights, prior.sd, y) {
    if (!is.null(etastart)) {
        eta <- etastart
    }
    else if (!is.null(start)) {
        start <- as.matrix (start)
        eta <- drop (x.nobs %*% start) + offset 
        ## TODO: should be able to drop this. coefold is the original start.
        #coefold <- start
    }
    else {
        eta <- family$linkfun(mustart)
    }
    if (family$family %in% c("poisson", "binomial")) {
        dispersion <- 1
    } else{
        dispersion <- var.y / 10000
    }
    
    mu <- family$linkinv(eta)
    mu.eta.val <- family$mu.eta(eta)
    dev <- sum (family$dev.resids (y, mu, weights))
    
    list (eta=eta,
            mu=mu,
            mu.eta.val=mu.eta.val,
            varmu=family$variance(mu),
            good=(weights > 0) & (mu.eta.val != 0),
            dispersion=dispersion,
            dev=dev,
            conv=FALSE,
            boundary=FALSE,
            prior.sd=prior.sd) 
}

.bayesglm.fit.loop.validateInputs <- function(etastart, start, nvars, xnames) {
    if (is.null(etastart) & !is.null(start) & length(start) != nvars) {
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                        nvars, paste(xnames, collapse = ", ")), domain = NA)
    }
}


.bayesglm.fit.loop.validateState <- function (state, family, control, iter, dispersionold, devold) {
    if (all(!state$good)) {
        warning("no observations informative at iteration ", state$iter)
        return (FALSE)
    }
    if (is.null(state$fit) == FALSE && any (!is.finite(state$fit$coefficients))) {
        warning("non-finite coefficients at iteration ", iter)
        return (FALSE)
    }
    if (any(is.na(state$varmu[state$good]))){
        stop("NAs in V(mu)")
    }
    if (any(state$varmu[state$good] == 0)){
        stop("0s in V(mu)")
    }
    
    if (any(is.na(state$mu.eta.val[state$good]))){
        stop("NAs in d(mu)/d(eta)")
    }
    
    if (iter > 1 & abs(state$dev - devold)/(0.1 + abs(state$dev)) <  control$epsilon & 
            abs(state$dispersion - dispersionold)/(0.1 + abs(state$dispersion)) < control$epsilon) {
        return (FALSE)
    }
    return (TRUE)
}

.bayesglm.fit.loop.updateState <- function (state, priors, family, #fortran.call.parameters, 
        offset, weights,
        y, x, x.nobs, nvars, nobs,
        intercept, scaled, control) {
    z <- (state$eta[state$good] - offset[state$good]) + (y[state$good] - state$mu[state$good]) / state$mu.eta.val[state$good]
    z.star <- c(z, priors$mean)
    
    w <- sqrt((weights[state$good] * state$mu.eta.val[state$good])^2/state$varmu[state$good])
    w.star <- c(w, sqrt(state$dispersion)/state$prior.sd)
    good.star <- c(state$good, rep(TRUE, nvars))
    fit <- lm.fit(x = as.matrix(x[good.star, ])*w.star, y = z.star*w.star)
 
    
    #fit <- .Fortran("dqrls", 
#            qr = x[good.star, ] * w.star, 
#            n = sum (good.star), 
#            p = nvars, 
#            y = w.star * z.star, 
#            ny = fortran.call.parameters$ny, 
#            tol = fortran.call.parameters$tol, 
#            coefficients = fortran.call.parameters$coefficients, 
#            residuals = fortran.call.parameters$residuals, 
#            effects = fortran.call.parameters$effects, 
#            rank = fortran.call.parameters$rank, 
#            pivot = fortran.call.parameters$pivot, 
#            qraux = fortran.call.parameters$qraux,
#            work = fortran.call.parameters$work, 
#            PACKAGE = "base")
    

        
    state$prior.sd <- priors$scale
    if (all (priors$df==Inf) == FALSE) {
      colMeans.x <- colMeans (x.nobs)  
      centered.coefs <- fit$coefficients
      if(NCOL(x.nobs)==1){
        V.coefs <- chol2inv(fit$qr$qr[1:nvars])
      }else{
        V.coefs <- chol2inv(fit$qr$qr[1:nvars, 1:nvars, drop = FALSE])
      }
      diag.V.coefs <- diag(V.coefs)
      sampling.var <- diag.V.coefs
            # DL: sampling.var[1] <- crossprod(colMeans.x, V.coefs) %*% colMeans.x
      if(intercept & scaled){      
        centered.coefs[1] <- sum(fit$coefficients*colMeans.x)
        sampling.var[1] <- crossprod (crossprod(V.coefs, colMeans.x), colMeans.x)
      }
      sd.tmp <- ((centered.coefs - priors$mean)^2 + sampling.var * state$dispersion + priors$df * priors$scale^2)/(1 + priors$df)
      sd <- sqrt(sd.tmp)
      state$prior.sd[priors$df != Inf] <- sd[priors$df != Inf]
    }
    
    if(NCOL(x.nobs)==1){
      predictions <- x.nobs * fit$coefficients
    }else{
      predictions <- x.nobs %*% fit$coefficients
    }
    
    if (!(family$family %in% c("poisson", "binomial"))) {
      if (exists ("V.coefs") == FALSE) {
        if(NCOL(x.nobs)==1){
          V.coefs <- chol2inv(fit$qr$qr[1:nvars])            
        }else{
          V.coefs <- chol2inv(fit$qr$qr[1:nvars, 1:nvars, drop = FALSE])
        }
      }
      #mse.resid <- mean((w * (z - x.nobs %*% fit$coefficients))^2) ## LOCAL VARIABLE
      #mse.resid <- mean ( (fit$y[1:nobs] - w * predictions)^2)
      ## mse.uncertainty <- mean(diag(x.nobs %*% V.coefs %*% t(x.nobs))) * state$dispersion
      mse.resid <- mean ( ((z.star*w.star)[1:sum(state$good)] - w * predictions[state$good,])^2)       
      mse.uncertainty <- max (0, mean(rowSums(( x.nobs %*% V.coefs ) * x.nobs)) * state$dispersion) #faster  ## LOCAL VARIABLE
      state$dispersion <- mse.resid + mse.uncertainty
    }
   
    state$eta <- drop(predictions)
    state$mu <- family$linkinv(state$eta + offset)
    state$mu.eta.val <- family$mu.eta(state$eta)
    
    dev <- sum (family$dev.resids (y, state$mu, weights))
 
    if (!is.finite (dev) || !isTRUE(family$valideta(state$eta) && family$validmu(state$mu))) {
        if (!is.finite(dev)) {
            warning("step size truncated due to divergence", call. = FALSE)
        } else if (!isTRUE(family$valideta(state$eta) && family$validmu(state$mu))) {
            warning("step size truncated: out of bounds", call. = FALSE)
        } 
         
        ii <- 1
        while (ii <= control$maxit & !is.finite (dev)) {
            ii <- ii + 1
            start <- (predictions + state$start) / 2
            state$mu <- family$linkinv(state$eta + offset)
            dev <- sum (family$dev.resids (y, state$mu, weights))
        }
        if (ii > control$maxit) {
            stop("inner loop 1; cannot correct step size")
        }
        
        ii <- 1
        while (ii <= control$maxit &  !isTRUE(family$valideta(state$eta) && family$validmu(state$mu))) {
            ii <- ii + 1
            start <- (predictions + state$start) / 2
            state$mu <- family$linkinv(state$eta + offset)
        }
        if (ii > control$maxit) {
            stop("inner loop 2; cannot correct step size")
        }
        
        
        state$boundary <- TRUE
        if (control$trace){
            cat("Step halved: new deviance =", dev, "\n")
        }
    }
    
    list (eta=state$eta,
            mu=state$mu,
            mu.eta.val=state$mu.eta.val,
            varmu=family$variance(state$mu),
            good=(weights > 0) & (state$mu.eta.val != 0),
            dispersion=state$dispersion,
            dev=dev,
            fit=fit,
            conv=FALSE,
            boundary=state$boundary,
            prior.sd=state$prior.sd,
            z=z,
            w=w) 
}


.bayesglm.fit.loop.initializePriors <- function (prior.mean, prior.mean.for.intercept, prior.scale, prior.df, prior.df.for.intercept) {
    list(mean=prior.mean,
            mean.for.intercept=prior.mean.for.intercept,
            scale=prior.scale,
            df=prior.df,
            df.for.intercept=prior.df.for.intercept)
}

.bayesglm.fit.loop.print <- function (state, priors, family, print.unnormalized.log.posterior, intercept, y) {
    if (print.unnormalized.log.posterior && family$family == "binomial") {
        logprior <- if(intercept) {
                    sum( dt( state$fit$coefficients[-1], priors$df , priors$mean,log = TRUE ) )
                    + dt( state$fit$coefficients[1], priors$df.for.intercept, priors$mean.for.intercept,
                            log = TRUE )
                }
                else {
                    sum( dt( state$fit$coefficients, priors$df , priors$mean, log = TRUE ) )
                }
        invlogit <- function (x) {
          1/(1+exp(-x))
        }
        #xb <- invlogit( x.nobs %*% coefs.hat )
        xb <- invlogit (state$eta - offset)
        loglikelihood <- sum( log( c( xb[ y == 1 ], 1 - xb[ y == 0 ] ) ) )
        cat( "log prior: ", logprior, ", log likelihood: ", loglikelihood, ",
                        unnormalized log posterior: ", loglikelihood +logprior, "\n" ,sep="")
    }
}


.bayesglm.fit.loop.createAuxillaryItems <- function (state, nvars, nobs, xnames, offset) {
    if (state$fit$rank < nvars) {
        state$fit$coefficients[seq(state$fit$rank + 1, nvars)] <- NA
    }
    residuals <- rep.int(NA, nobs)
    residuals[state$good] <- state$z[state$good] - (state$eta[state$good] - offset[state$good])
    
    state$fit$qr$qr <- as.matrix(state$fit$qr$qr)
    
    nr <- min(sum(state$good), nvars)
    if (nr < nvars) {
      Rmat <- diag(x = 0, nvars)
      Rmat[1:nr, 1:nvars] <- state$fit$qr$qr[1:nr, 1:nvars, drop=FALSE]
    }
    else{
      if(NCOL(state$fit$qr$qr)==1){
        Rmat <- state$fit$qr$qr[1:nvars]
      }else{
        Rmat <- state$fit$qr$qr[1:nvars, 1:nvars, drop=FALSE]
      }
    }
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names (state$fit$coefficients) <- xnames
    colnames(state$fit$qr$qr) <- xnames
    dimnames(Rmat) <- list(xnames, xnames)
    
    list (residuals=residuals,
            Rmat=Rmat,
            state=state)    
}

.bayesglm.fit.loop.printWarnings <- function (Warning, state, family) {
    if(Warning){
        if (!state$conv){
            warning("algorithm did not converge")
        }
        if (state$boundary){
            warning("algorithm stopped at boundary value")
        }
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(state$mu > 1 - eps) || any(state$mu < eps)) {
                warning("fitted probabilities numerically 0 or 1 occurred")
            }
        }
        if (family$family == "poisson") {
            if (any(state$mu < eps)){
                warning("fitted rates numerically 0 occurred")
            }
        }
    }
}

.bayesglm.fit.loop.main.ideal <- function (control, x, y, nvars, nobs, weights, offset,
        intercept, scaled, 
        start, etastart, mustart,
        family,
        prior.mean, prior.mean.for.intercept, prior.scale, prior.df, prior.df.for.intercept,
        print.unnormalized.log.posterior,
        Warning) {
    xnames <- dimnames(x)[[2]]
    x.nobs <- x[1:nobs, ,drop=FALSE]
    
    .bayesglm.fit.loop.validateInputs (etastart, start, nvars, dimnames(x)[[2]])
    
    # invalid starting point (should be moved out of here):
    # is.null(start)==false & length(start) != nvars
    
    priors <- .bayesglm.fit.loop.initializePriors (prior.mean, prior.mean.for.intercept, prior.scale, prior.df, prior.df.for.intercept)    
    
    state <- .bayesglm.fit.loop.initializeState (start, etastart, mustart, offset, x.nobs, var(y), nvars, family, weights, priors$scale, y)
    #core elements of state:
    # good: derived from eta
    # eta
    # x, y, nvars, nobs: invariants
    
    #fortran.call.parameters <- .bayesglm.fit.loop.initialMemoryAllocation (control$epsilon, nvars, sum (state$good))
    
    for (iter in 1:control$maxit) {
        dispersionold <- state$dispersion
        devold <- state$dev
        state <- .bayesglm.fit.loop.updateState (state, priors, family, #fortran.call.parameters, 
                offset, weights,
                y, x, x.nobs, nvars, nobs,
                intercept, scaled, control) 
        if (control$trace){
            cat("Deviance =", state$dev, "Iterations -", iter, "\n")
        }
        
        if (.bayesglm.fit.loop.validateState(state, family, control, iter, dispersionold, devold) == FALSE) {
            state$conv <- TRUE
            break
        }
        .bayesglm.fit.loop.print(state, priors, family, print.unnormalized.log.posterior, intercept, y)
    }
    
    .bayesglm.fit.loop.printWarnings(Warning, state, family)
    output <- .bayesglm.fit.loop.createAuxillaryItems(state, nvars, nobs, xnames, offset)
    ## residuals=residuals
    ## Rmat=Rmat
    ## state
    
    list (fit=output$state$fit,
            good=output$state$good,
            z=output$state$z,
            ## z.star=output$state$z.star,
            w=output$state$w,
            ## w.star=output$state$w.star,
            ngoodobs=sum (output$state$good),
            prior.scale=priors$scale,
            prior.sd=output$state$prior.sd,
            eta=output$state$eta,
            mu=output$state$mu,
            dev=output$state$dev,
            dispersion=output$state$dispersion,
            dispersionold=dispersionold,
            start=output$state$fit$coefficients,
            coef=output$state$fit$coefficients,
            devold=devold,
            conv=output$state$conv,
            iter=iter,
            boundary=output$state$boundary,
            Rmat=output$Rmat,
            residuals=output$residuals)
}


.bayesglm.fit.cleanup <- function (ynames, residuals, mu, eta, nobs, wt, w, good, weights, linkinv, dev.resids, y, intercept, fit, offset, EMPTY, n, dev, aic) {
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    wtdmu <- if (intercept){
                sum(weights * y)/sum(weights)
            }
            else{
                linkinv(offset)
            }
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY){
                0
            }
            else{
                fit$rank
            }
    resdf <- n.ok - rank
    resdf <- n.ok
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    
    list (residuals=residuals,
            mu=mu,
            eta=eta,
            wt=wt,
            weights=weights,
            y=y,
            wtdmu=wtdmu,
            nulldev=nulldev,
            n.ok=n.ok,
            nulldf=nulldf,
            rank=rank,
            resdf=resdf,
            aic.model=aic.model)
}

getQr <- function(x, ...){
  if (is.null(r <- x$qr)) 
        stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
  r
}


## from model.matrix.bayes.R in arm

#setMethod("model.matrix.bayes", signature(object = "bayesglm"),
model.matrixBayes.internal <- function(object, data = environment(object),
                              contrasts.arg = NULL, xlev = NULL, keep.order=FALSE, drop.baseline=FALSE,...)
{
  #class(object) <- c("terms", "formula")
  t <- if( missing( data ) ) { 
    terms( object ) 
  }else{ 
    terms.formula(object, data = data, keep.order=keep.order) 
  }
  attr(t, "intercept") <- attr(object, "intercept")
  if (is.null(attr(data, "terms"))){ 
    data <- model.frame(object, data, xlev=xlev) 
  }else {
    reorder <- match(sapply(attr(t,"variables"), deparse, width.cutoff=500)[-1], names(data))
    if (any(is.na(reorder))) {
      stop( "model frame and formula mismatch in model.matrix()" ) 
    }
    if(!identical(reorder, seq_len(ncol(data)))) {
      data <- data[,reorder, drop = FALSE] 
    }
  }
  int <- attr(t, "response")
  if(length(data)) {      # otherwise no rhs terms, so skip all this
    
    if (drop.baseline){
      contr.funs <- as.character(getOption("contrasts"))
    }else{
      contr.funs <- as.character(list("contr.bayes.unordered", "contr.bayes.ordered"))
    }
    
    namD <- names(data)
    ## turn any character columns into factors
    for(i in namD)
      if(is.character( data[[i]] ) ) {
        data[[i]] <- factor(data[[i]])
        warning( gettextf( "variable '%s' converted to a factor", i ), domain = NA)
      }
    isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)        
    isF[int] <- FALSE
    isOF <- vapply(data, is.ordered, NA)
    for( nn in namD[isF] )            # drop response
      if( is.null( attr( data[[nn]], "contrasts" ) ) ) {
        contrasts( data[[nn]] ) <- contr.funs[1 + isOF[nn]]
      }
    ## it might be safer to have numerical contrasts:
    ##    get(contr.funs[1 + isOF[nn]])(nlevels(data[[nn]]))
    if ( !is.null( contrasts.arg ) && is.list( contrasts.arg ) ) {
      if ( is.null( namC <- names( contrasts.arg ) ) ) {
        stop( "invalid 'contrasts.arg' argument" )
      }
      for (nn in namC) {
        if ( is.na( ni <- match( nn, namD ) ) ) {
          warning( gettextf( "variable '%s' is absent, its contrast will be ignored", nn ), domain = NA )
        }
        else {
          ca <- contrasts.arg[[nn]]
          if( is.matrix( ca ) ) {
            contrasts( data[[ni]], ncol( ca ) ) <- ca
          }
          else { 
            contrasts( data[[ni]] ) <- contrasts.arg[[nn]]
          }
        }
      }
    }
  } else {               # internal model.matrix needs some variable
    isF  <-  FALSE
    data <- data.frame(x=rep(0, nrow(data)))
  }
  #ans  <- .Internal( model.matrix( t, data ) )
  ans  <- model.matrix.default(object=t, data=data)
  cons <- if(any(isF)){
    lapply( data[isF], function(x) attr( x,  "contrasts") ) 
  }else { NULL }
  attr(ans, "contrasts" ) <- cons
  ans
}

