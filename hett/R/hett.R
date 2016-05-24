"tlm" <-
function(lform = formula(data), sform =  ~ 1, data = sys.parent(), subset = NULL, contrasts = NULL, na.action = na.fail, start = NULL, control = tlm.control(...), obs = FALSE, estDof = FALSE, ...)
{
  ## A function for modelling the location and scale parameters of the Students t..
  ## This will only handle a specific case where conditional on \omega, y = Xb + e as the location model
  ## where e \sim N(0, \sigma^2/\omega) and log(\sigma^2) = Z\lambda + log(\omega), where \omega
  ## is distributed as Gamma(\nu/2,\nu/2). Marginally, y \sim t(Xb, \sigma^2, \nu). Reference (Lange, Little and Taylor, 1989)) 
  ##
  ## Arguments :
  ##       lform   : location formula
  ##       sform   : scale formula
  ##       start   : starting values (list of beta,lambda,dof,omega)
  ##       control : control prameters for iterative scheme
  ##       obs     : TRUE, FALSE. Use the expected information for estimation or reweighted least squares estimation
  ##       estDof  : TRUE, FALSE. Estimate or fix the degrees of freedom parameter

  ## Comments :
  ## 1) Fitting a model with constant scale requires no argument for "sform=".
  ## 2) For the initial dof in the start argument
  ##               : a) If estimating the degrees of freedom then this value is the intial (default = 3) only.
  ##                 b) This is the value of the degrees of freedom if it is not being estimated (default = 3). 
  ## 2) To obtain the correct standard errors for the estimates use "obs = FALSE" so
  ##    the fit uses the expected information for its iteration of the location.
  ## 3) At the moment testing between models should be done using difference in
  ##    likelihoods or the score test rather than using the Pr(z-value) from summary().
  ##    This can be done using the "tscore" function.
  
  tcall <- match.call()	##
  ## Setting up the location and scale model and
  ## extracting the appropriate information.
  ##
  sform[3] <- sform[2]
  sform[2] <- lform[2]
  fixedFrameArgs <- list(lform, na.action = na.action, contrasts = contrasts, data = data)
  dispFrameArgs <- list(sform, na.action = na.action, data = data)
  if(!missing(subset)) {
    fixedFrameArgs[["subset"]] <- subset
    dsipFrameArgs[["subset"]] <- subset
  }
  fixedFrame <- do.call("model.frame", fixedFrameArgs)
  dispFrame <- do.call("model.frame", dispFrameArgs)
  X <- model.matrix(lform, fixedFrame)
  y <- model.extract(fixedFrame, response)
  Z <- model.matrix(sform, dispFrame)
  if(!is.null(start)) {
    if(is.null(namStart <- names(start)))
      stop("Starting values must have names for assignation")
    namList <- c("beta", "lambda", "dof", "omega")
    if(!all(pmatch(namStart, namList, nomatch = 0)))
      stop("Names of starting values must be \"beta\", \"lambda\", \"dof\" or \"omega\"")
  }
  else namStart <- ""
  if((bets <- pmatch("beta", namStart, nomatch = 0)) != 0) {
    fixedCoef <- start[[bets]]
    sqResid <- (y - crossprod(t(X), fixedCoef))^2
  }
  else {
    initFit <- lm.fit(X, y, method = "qr")
    fixedCoef <- coef(initFit)
    sqResid <- resid(initFit)^2
  }
  dispFamily <- Gamma(link = "log") 
  if((lams <- pmatch("lambda", namStart, nomatch = 0)) != 0) {
    dispCoef <- start[[lams]]
    sigmaI <- dispFamily$linkinv(crossprod(t(Z), dispCoef))
  }
  else {
    dispFit <- lm.fit(Z, dispFamily$linkfun((sqResid + (sqResid == 0)/6)))
    zLambda <- dispFit$fitted.values
    dispCoef <- coef(dispFit)
    sigmaI <- dispFamily$linkinv(zLambda)
  }
  if((dofs <- pmatch("dof", namStart, nomatch = 0)) != 0){
    dof <- start[[dofs]]
    if(length(dof) != 1) 
      stop(paste("initial value for degrees of freedom must be of length 1"))
  }
  else   
    dof <- 3
  if((oms <- pmatch("omega", namStart, nomatch = 0)) != 0){
    randCoef <- start[[oms]]
    if(length(randCoef) != length(y))
      stop(paste("Starting values for scale random effects must be the same length as y"))
  }
  else 
    randCoef <- (dof + 1)/(dof + sqResid/sigmaI)
  n <- length(y)
  ## control parameters for algorithm
  epsilon <- control$epsilon
  maxit <- control$maxit
  trace <- control$trace
  verbose <- control$verboseLev
  iter <- 0
  initTime <- proc.time()[1]
  ##
  ## iterative scheme in 4 parts : location -> fixed scale -> random scale -> degrees of freedom
  ##
  max.fit <- expression({
    const <-  - n * lgamma(1/2) - n * lgamma(dof/2) - (n/2) * log(dof) + n * lgamma((dof + 1)/2)
    logLik <- const - (sum(log(sigmaI)))/2 - sum(((dof + 1)/2) * (log(1 + (sqResid/sigmaI)/dof)))
    logLikTemp <- logLik + 1
    orthoComp <- apply(crossprod(solve(crossprod(Z)), t(Z)),1,sum)
    orthoCoef <- dispCoef - orthoComp 
    while(abs(logLik - logLikTemp) > epsilon && iter < maxit) {    
      ##            
      ## location parameter likelihood
      ##
      if(!obs) {
        mscale <- (dof + 1)/(dof + 3)
        meanY <-  (y - X %*% fixedCoef)*(randCoef/mscale) + crossprod(t(X), fixedCoef)
        meanFit <- lm.wfit(X, meanY, rep(mscale, length(y))/sigmaI, method = "qr")
      }
      else meanFit <- lm.wfit(X, y, randCoef/sigmaI, method = "qr")
      fitm <- meanFit$fitted.values
      fixedCoef <- coef(meanFit)
      sqResid <- (as.vector(y - crossprod(t(X), fixedCoef)))^2
      if(trace) cat("\nLocation parameters :", 
                    format(round(fixedCoef, 4)))	
      ##          
      ## scale parameter estimation
      ##
      tscale <- (1 - 3/(dof + 3))
      scale <- (1/tscale)
      dispY <- scale * (((sqResid * randCoef)/sigmaI) - rep(1, length(y))) + crossprod(t(Z), dispCoef)
      dispFit <- lm.wfit(Z, dispY, rep(tscale, length(y)), method = "qr")
      dispCoef <- coef(dispFit)
      sigmaI <- dispFamily$linkinv(dispFit$fitted.values)
      orthoI <- as.vector(sigmaI/exp(Z %*% (orthoComp*(2*(log(dof) - log(dof + 1))))))
      if(trace) cat("\nScale parameters :", format(round(
                                                              dispCoef, 4)))
      ##
      ## random component maximum likelihood 
      ##
      randCoef <- (dof + 1)/(dof + sqResid/sigmaI)
      if(trace && verbose == 2)
        cat("\nRandom scale effects :", format(round(randCoef, 4)))
      ##
      ## degrees of freedom estimation
      ##
      if(estDof) {
        opt <- nlm(dof.profile, dof, n = n, sqResid = sqResid, orthoI = orthoI, X = X, Z = Z)
        dof <- opt$estimate
        if(trace)
          cat("\ndegrees of freedom is ", dof)
        ## expected information for degrees of freedom standard error. 
        expinfDof <-   n*trigamma(dof/2)/4 - n*trigamma((dof + 1)/2)/4 - n*(dof + 5)/(2 * dof * (dof + 1)*(dof + 3)) - (2/(dof*(dof + 3)*(dof + 1)^2))*(t(rep(1, length(y))) %*% Z %*% orthoComp)
      }
      logLikTemp <- logLik
      if(estDof) {
        logLik <- - opt$minimum
        dofse <- 1/sqrt(expinfDof)
      }
      else {
        logLik <- const - sum(log(sigmaI))/2 - sum(((dof + 1)/2) * (log(1 + (sqResid/sigmaI)/dof)))
        dofse <- NA
      }
      iter <- iter + 1
      if(trace)
        cat("\nHeteroscedastic t: - 2 x Maximum Likelihood iteration", iter, " : ", format(round(-2 * logLik, 4)), "\n")
    }})
  estMethod <-  eval(max.fit)        
  ##
  ## Coverged or not?
  ##                 
  endTime <- proc.time()[1] - initTime
  if(trace && iter >= maxit) cat("\nMaximum", maxit, 
                "iterations reached in", endTime, "\n")	                                      
  meanFit$method <- "maximum likelihood"
  random <- randCoef
  meanFit$formula <- lform
  meanFit$terms <- attr(fixedFrame, "terms")
  meanFit$iter <- iter
  meanFit$call <- tcall
  meanFit$residuals <- as.vector(y - crossprod(t(X), fixedCoef)) 
  dispFit$formula <- sform
  dispFit$terms <- attr(dispFrame, "terms")
  dispFit$iter <- iter
  dispFit$fitted.values <- sigmaI
  dispFit$call <- tcall
  res <- list(loc.fit = meanFit, scale.fit = dispFit, random = 
              random, dof = dof, dofse = dofse, iter = iter, logLik = 
              logLik, endTime = endTime)
  class(res) <- c("tlm", "glm", "lm")
  res
}

"summary.tlm" <-
function(object, correlation = FALSE, ...)
{
  locOut <- object$loc.fit
  scaleOut <- object$scale.fit
  locSum <- tsum(locOut, dispersion = 1, correlation = correlation, ...)
  scaleSum <- tsum(scaleOut, dispersion = 2, correlation = correlation, ...)
  structure(list(loc.summary = locSum, scale.summary = scaleSum, iter = 
                 object$iter, dof = object$dof, dofse = object$dofse,
                 logLik = object$logLik, method = locOut$method, 
                 estTime = object$endTime), class = c("summary.tlm"))
}

"print.tlm" <-
function(x, ...)
  {
    cat("\nCall: ")
    cat(format(x$loc.fit$call), "\n\n")
    cat("Location Coefficients:\n")
    print.default(format(x$loc.fit$coefficients, digits = options()$digits), print.gap = 2, 
                  quote = FALSE)
    cat("Scale Coefficients:\n")
    print.default(format(x$scale.fit$coefficients, digits = options()$digits), print.gap = 2, 
                  quote = FALSE)
    cat("Degrees of freedom:\n")
    cat(format(x$dof), "\n")
    cat("Log-Likelihood\n")
    cat(format(x$logLik), "\n")    
  }

"print.summary.tlm" <-
function(x, ...)
{
  cat("Location model :\n")
  print(x$loc.summary, scale = FALSE, ...)
  cat("\nScale Model :\n")
  print(x$scale.summary, scale = TRUE, ...)
  if(is.null(x$dofse))
    cat("\nFixed degrees of freedom parameter: ", format(x$dof))
  else {
    cat("\nEst. degrees of freedom parameter: ", format(x$dof))
    cat("\nStandard error for d.o.f: ", format(x$dofse))
  }
  cat("\nNo. of iterations of model :", format(x$iter), "in", format(
                                                                          round(x$estTime, options()$digits)))
  cat("\nHeteroscedastic t Likelihood :", format(x$logLik), "\n")
}

"tsum" <-
function (object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, 
            ...)
  ## Function edited from summary.glm()
  ## Sets up the summary information for the location and scale model for the t
  ## Produces the correct standard errors and correlations between coefficients.
{
  Qr <- object$qr
  est.disp <- FALSE
  df.r <- object$df.residual
  if (dispersion == 1)
    est.disp <- TRUE
  p <- object$rank
  p1 <- 1:p
  coef.p <- object$coefficients[Qr$pivot[p1]]
  covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
  covmat <- dispersion * covmat.unscaled
  var.cf <- diag(covmat)
  s.err <- sqrt(var.cf)
  tvalue <- coef.p/s.err
  dn <- c("Estimate", "Std. Error")
  if (!est.disp) {
    pvalue <- 2 * pnorm(-abs(tvalue))
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", 
                                                  "Pr(>|z|)"))
  }
  else if (df.r > 0) {
    pvalue <- 2 * pt(-abs(tvalue), df.r)
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", 
                                                  "Pr(>|t|)"))
  }
  else {
    coef.table <- cbind(coef.p, Inf)
    dimnames(coef.table) <- list(names(coef.p), dn)
  }
  ans <- c(object[c("call", "terms", "df.residual", "iter")], list(resid = residuals(object), coefficients = coef.table, dispersion = dispersion, df = c(object$rank, df.r), cov.unscaled = covmat.unscaled, cov.scaled = covmat))
  if (correlation) {
    dd <- sqrt(diag(covmat.unscaled))
    ans$correlation <- covmat.unscaled/outer(dd, dd)
    ans$symbolic.cor <- symbolic.cor
  }
  class(ans) <- "tsum"
  return(ans)
}

"print.tsum" <-
function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
            signif.stars = getOption("show.signif.stars"), scale = TRUE, ...) 
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Residuals: \n")
  if (x$df.residual > 5) {
    x$resid <- quantile(x$resid, na.rm = TRUE)
    names(x$resid) <- c("Min", "1Q", "Median", "3Q", 
                        "Max")
  }
  print.default(x$resid, digits = digits, na = "", 
                print.gap = 2)
  cat("\nCoefficients:\n")
  printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, 
                ...)
  if(scale)
    cat("\n(Scale parameter taken to be ", 
        format(x$dispersion),")\n")
  else
    cat("\n(Scale parameter(s) as estimated below)\n")  
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.col = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2, 
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}

"tlm.control" <-
function(epsilon = 1e-07, maxit = 50, trace = FALSE, verboseLev = 1)
{
  ## tdispersion control parameters. 
  if(epsilon <= 0) {
    warning("the value of epsilon supplied is zero or negative;\nthe default value of 1e-7 was used instead"
            )
    epsilon <- 1e-07
  }
  if(maxit < 1) {
    warning("the value of maxit supplied is zero or negative;\nthe default value of 50 was used instead"
            )
    maxit <- 50
  }
  if(verboseLev < 0 || verboseLev > 2)
    warning("the value of the verbose level is not in the expected range. Use \"verboseLev = 1\" for simple output and \"verboseLev = 2\" for advanced output"
            )
  list(epsilon = epsilon, maxit = maxit, trace = trace, 
       verboseLev = verboseLev)
}

"dof.profile" <-
function(dof, n, sqResid, orthoI, X, Z)
{
  n <-  length(sqResid)
  sigmaI <- as.vector(orthoI*exp(Z %*% crossprod(solve(crossprod(Z)), t(Z)) %*% rep(2*(log(dof) - log(dof + 1)), n)))
  - (- n * lgamma(1/2) - n * lgamma(dof/2) - (n/2) * log(dof) + n * lgamma((dof + 1)/2) -
     (sum(log(sigmaI)))/2 - sum(((dof + 1)/2) * (log(1 + (sqResid/sigmaI)/dof))))
}

"tscore" <-
function(..., data = NULL, scale = FALSE)
{
  ## General score test for the t-dispersion model for both scale and location
  listm <- list(...)
  if(is.null(data))
    stop("Score function requires the data to reproduce some model elements")
  assign("data", data, envir = .GlobalEnv)
  dofse <- unlist(lapply(listm, function(el) el$dofse))
  if(any(is.na(dofse))){
    if(!all(is.na(dofse)))
      stop("All models must have estimated or fixed degrees of freedom\n")
  }
  fullm <- do.call("model.frame", list(listm[[length(listm)]]$loc.fit$formula, data = data))
  y <- model.extract(fullm, response)  
  k <- length(listm)
  ss <- pv <- df <- vector()
  if(scale)
    Xl <- lapply(listm, function(el, data) model.matrix(el$scale.fit$formula,  data = data), data = data)
  else
    Xl <- lapply(listm, function(el, data) model.matrix(el$loc.fit$formula, data = data), data = data)
  k <- length(listm)
  nams <- vector()
  j <- 1
  while(k > 1){
    X2 <- as.matrix(Xl[[k]][, !pmatch(dimnames(Xl[[k]])[[2]], dimnames(Xl[[k - 1]])[[2]], nomatch = 0)])
    random <- listm[[k - 1]]$random
    ssig <- (listm[[k - 1]]$scale.fit$fitted.values)^(1/2)
    Om <- random/(ssig^2)
    if(scale){
      Q <- qr.Q(qr(Xl[[k - 1]]))
      P <- - crossprod(t(Q))
      diag(P) <- 1 + diag(P)
      ss[j] <- ((listm[[k - 1]]$dof + 3)/(2 * listm[[k - 1]]$dof))*(t((resid(listm[[k - 1]]$loc.fit)^2)*Om - 1) %*% X2 %*% solve(t(X2) %*% P %*% X2)
                                                                    %*% t(X2) %*% (Om*(resid(listm[[k - 1]]$loc.fit)^2) - 1))
      df[j] <- listm[[k - 1]]$scale.fit$df.residual - listm[[k]]$scale.fit$df.residual
      nams[j] <- paste(paste(deparse(listm[[k - 1]]$scale.fit$formula[[3]]), "vs", sep = " "), deparse(listm[[k]]$scale.fit$formula[[3]]), sep = " ")
    }
    else { 
      Q <- qr.Q(qr(ssig*Xl[[k - 1]]))
      P <- diag(1/ssig^2) - diag(1/ssig) %*% crossprod(t(Q)) %*% diag(1/ssig) 
      ss[j] <-  ((listm[[k - 1]]$dof + 3)/(listm[[k - 1]]$dof + 1))*(t(resid(listm[[k - 1]]$loc.fit)*Om) %*% X2 %*% solve(t(X2) %*% P %*% X2)
                                                                     %*% t(X2) %*% (Om*resid(listm[[k - 1]]$loc.fit)))
      df[j] <- listm[[k - 1]]$loc.fit$df.residual - listm[[k]]$loc.fit$df.residual
      nams[j] <- paste(paste(deparse(listm[[k - 1]]$loc.fit$formula[[3]]), "vs", sep = " "), deparse(listm[[k]]$loc.fit$formula[[3]]), sep = " ")
    }
    pv[j] <- 1 - pchisq(ss[j], df[j])    
    k <- k - 1
    j <- j + 1
  }
  tab <- cbind(df, ss, pv)
  colnames(tab) <- c("df", "Score Stat.", "P-value")
  rownames(tab) <- nams
  cat("Heteroscedastic t Regression:\n")
  if(scale)
    cat("\nScale score statistic tests:\n\n")
  else
    cat("\nLocation score statistic tests:\n\n")
  printCoefmat(tab)
}

