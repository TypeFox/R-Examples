onls <- function(formula, data = parent.frame(), start, jac = NULL, 
                 control = nls.lm.control(), lower = NULL, 
                 upper = NULL, trace = FALSE, subset, weights, na.action, 
                 window = 12, extend = c(0.2, 0.2), fixed = NULL, 
                 verbose = TRUE, ...) 
{
  options(warn = -1)
  formula <- as.formula(formula)
  if (!is.list(data) && !is.environment(data)) 
    stop("'data' must be a list or an environment")
  mf <- match.call()
  if (!is.null(fixed) & length(fixed) != length(start))
    stop("'fixed' must have same length as 'start'!")
  
  ## remove 'extend', 'fixed' , 'window' and 'verbose' from call, so 'model.frame' throws no error
  mf <- as.list(mf)
  EXTEND <- mf$extend  
  mf$extend <- NULL
  FIXED <- fixed
  mf$fixed <- NULL
  mf$verbose <- NULL
  mf$window <- NULL
  mf <- as.call(mf)   
    
  varNames <- all.vars(formula)
  if (length(formula) == 2L) {
    formula[[3L]] <- formula[[2L]]
    formula[[2L]] <- 0
  }
  form2 <- formula
  form2[[2L]] <- 0
  
  ## get name of all RHS variables
  varNamesRHS <- all.vars(form2)
    
  ## get name of response variable
  varMatch <- match(varNamesRHS, varNames) 
  varNamesLHS <- varNames[-varMatch]
  
  mWeights <- missing(weights)
  if (trace) control$nprint <- 1
  pnames <- if (missing(start)) {
    if (!is.null(attr(data, "parameters"))) {
      names(attr(data, "parameters"))
    }
    else {
      cll <- formula[[length(formula)]]
      func <- get(as.character(cll[[1L]]))
      if (!is.null(pn <- attr(func, "pnames"))) 
        as.character(as.list(match.call(func, call = cll))[-1L][pn])
    }
  }
  else names(start)
    
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  if (length(pnames)) varNames <- varNames[is.na(match(varNames, pnames))]
  lenVar <- function(var) tryCatch(length(eval(as.name(var), data, env)), error = function(e) -1)
  
  if (length(varNames)) {
    n <- sapply(varNames, lenVar)
    if (any(not.there <- n == -1)) {      
      nnn <- names(n[not.there])      
      if (missing(start)) {
        warning("No starting values specified for some parameters.\n", 
                "Initializing ", paste(sQuote(nnn), collapse = ", "), 
                " to '1.'.\n", "Consider specifying 'start' or using a selfStart model")
        start <- as.list(rep(1, length(nnn)))        
        names(start) <- nnn
        varNames <- varNames[i <- is.na(match(varNames, nnn))]
        n <- n[i]
      }
      else stop("parameters without starting value in 'data': ", 
                paste(nnn, collapse = ", "))
    }
  }
  else {        
    if (length(pnames) && any((np <- sapply(pnames, lenVar)) == -1)) {
      message("fitting parameters ", paste(sQuote(pnames[np == -1]), 
                                           collapse = ", "), " without any variables")
      n <- integer()
    }
    else stop("no parameters to fit")
  }
    
  respLength <- length(eval(formula[[2L]], data, env))
    
  if (length(n) > 0L) {
    varIndex <- n%%respLength == 0
    
    if (is.list(data) && diff(range(n[names(n) %in% names(data)])) > 0) {
      mf <- data 
      if (!missing(subset)) warning("argument 'subset' will be ignored")
      if (!missing(na.action)) warning("argument 'na.action' will be ignored")
      if (missing(start)) start <- getInitial(formula, mf)
            
      startEnv <- new.env(hash = FALSE, parent = environment(formula))
      
      for (i in names(start)) assign(i, start[[i]], envir = startEnv)
      rhs <- eval(formula[[3L]], data, startEnv)
      n <- NROW(rhs)
      wts <- if (mWeights) rep(1, n)
      else eval(substitute(weights), data, environment(formula))
    }
    else {
      mf$formula <- as.formula(paste("~", paste(varNames[varIndex], 
                               collapse = "+")), env = environment(formula))
      mf$start <- mf$control <- mf$algorithm <- mf$trace <- mf$model <- NULL
      mf$lower <- mf$upper <- NULL
      mf[[1L]] <- as.name("model.frame")
      mf <- eval.parent(mf)
      n <- nrow(mf)
      mf <- as.list(mf)
      wts <- if (!mWeights) model.weights(mf)
      else rep(1, n)      
    }
    
    if (any(wts < 0 | is.na(wts))) stop("missing or negative weights not allowed")
  }
  else {
    varIndex <- logical()
    mf <- list(0)
    wts <- numeric()
  }
     
  if (missing(start)) start <- getInitial(formula, mf)
  for (var in varNames[!varIndex]) mf[[var]] <- eval(as.name(var), data, env)
  
  ## get name of predictor variable
  varNamesRHS <- varNamesRHS[varNamesRHS %in% varNames[varIndex]]
  
  mf <- c(mf, start)
  lhs <- eval(formula[[2L]], envir = mf)
  m <- match(names(start), names(mf))
  .swts <- if (!missing(wts) && length(wts)) sqrt(wts)
  
  ## set control arguments for 'nls.lm'
  control$maxiter <- 1024  
  control$maxfev <- 100000
      
  ## get x, y values and set border to optimize in,
  ## create new environment for storing X0 values
  PRED <- mf[[varNamesRHS]] # x   
  RESP <- mf[[varNamesLHS]] # y  
  LOWER <- min(PRED, na.rm = TRUE) - extend[1] * diff(range(PRED, na.rm = TRUE))
  UPPER <- max(PRED, na.rm = TRUE) + extend[2] * diff(range(PRED, na.rm = TRUE))
  envX0s <- new.env()    
  
  ## order predictor and response values
  ORD <- order(PRED)
  PRED <- PRED[ORD]
  RESP <- RESP[ORD]  
  
  ## usual NLS from nlsLM, objective function minimizes
  ## vertical sum-of-squares  
  nlsFCT <- function(par) {
    mf[m] <- par
    rhs <- eval(formula[[3L]], envir = mf)
    res <- lhs - rhs
    res <- .swts * res
    res
  }
  
  ## ordinary NLS, get fitted and residual values
  if (verbose) cat("Obtaining starting parameters from ordinary NLS...\n")
  
  NLS <- nls.lm(par = start, fn = nlsFCT, jac = NULL, control = control, 
                lower = lower, upper = upper, ...)
  if (verbose) {
    if (NLS$info %in% 1:4) cat("  Passed...\n", NLS$message, "\n")
    else cat("  Failed...", NLS$message, "\n")  
  }

  residNLS <- NLS$fvec
  fittedNLS <- RESP - residNLS
    
  ## orthogonal LS objective function minimizes
  ## euclidean distance => 'optimize'
  optFCT <- function(x, x2, y2, mf) {    
    mf2 <- mf   
    mf2[[varNamesRHS]] <- x    
    y <- eval(formula[[3L]], envir = mf2)    
    DIST <- sqrt((x - x2)^2 + (y - y2)^2)  
    DIST    
  } 
   
  ## pre-allocate vectors
  resid_o <- X0s <- numeric(length(RESP)) 
  LEN <- length(resid_o)
  
  ## convert start to numeric and define fixed parameters
  ## change starting value for 'ONLS'
  START <- unlist(start)
  SEL <- which(FIXED == TRUE)  
  NLS$par[SEL] <- START[SEL]
  
  ## orthogonal LS function, outer loop => LM,
  ## inner loop => 'optimize' for each x/y pair
  onlsFCT <- function(par, mf) {
    par[SEL] <- START[SEL]    
    mf[m] <- par    
        
    for (i in 1:LEN) {      
      ## define default intervals
      if (LEN <= 25) {
        INTERVAL <- c(LOWER, UPPER)
      } else {
        if (i <= window) INTERVAL <- c(LOWER, PRED[i + (window - 1)])
        else if (i >= LEN - (window - 1)) INTERVAL <- c(PRED[i - (window - 1)], UPPER)
        else INTERVAL <- c(PRED[i - (window - 1)], PRED[i + (window - 1)])        
      }         
            
      OPT <- optimize(optFCT, interval = INTERVAL, x2 = PRED[i], y2 = RESP[i], mf = mf,
                      tol = .Machine$double.eps^0.5) 
      resid_o[i] <- OPT$objective  
      X0s[i] <<- OPT$minimum                
    }
      
    res <- .swts * resid_o  
    res
  }
  
  ## do orthogonal fitting with starting parameters from NLS
  if (verbose) cat("Optimizing orthogonal NLS...\n")  
  ONLS <- nls.lm(par = NLS$par, fn = onlsFCT, jac = NULL, control = control, 
                 lower = lower, upper = upper, mf = mf)
    
  if (verbose) {
    if (ONLS$info %in% 1:4) cat("  Passed...", NLS$message, "\n\n")
    else cat("  Failed...", NLS$message, "\n") 
  }
  
  ## fitted, vertical and orthogonal residuals values 
  mf[m] <- ONLS$par
  fittedONLS <- eval(formula[[3L]], envir = mf)
  residONLS <- RESP - fittedONLS
  resid.o <- ONLS$fvec
  
  ## y0 values from orthogonal minimization
  mf4 <- mf  
  mf4[[varNamesRHS]] <- X0s
  Y0s <- eval(formula[[3L]], envir = mf4) 
  
  ## numeric gradient
  thisEnv <- environment()
  env <- new.env(hash = TRUE, parent = environment(formula))
  for (i in names(data)) assign(i, data[[i]], envir = env)
  PARS <- as.list(ONLS$par)  # converged parameters
  for (i in names(PARS)) {
    temp <- PARS[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp, envir = env)    
  }
  if (is.null(upper)) GRAD <- numericDeriv(formula[[3L]], pnames, env)
  else GRAD <- numericDeriv(formula[[3L]], pnames, env, ifelse(PARS < upper, 1, -1))
  QR <- qr(.swts * attr(GRAD, "gradient"))
  qrDim <- min(dim(QR$qr))
  if (QR$rank < qrDim) warning("Singular gradient matrix at initial parameter estimates")
  
  ## create output  
  if (ONLS$info %in% 1:4) isConv <- TRUE else isConv <- FALSE
  finIter <- ONLS$niter
  finTol <- control$ftol
  convInfo <- list(isConv = isConv, finIter = ONLS$niter, finTol = finTol, 
                   stopCode = ONLS$info, stopMessage = ONLS$message)
  OUT <- list(data = substitute(data), call = match.call(),
              convInfo = convInfo)  
  OUT$call$control <- control  
  OUT$call$trace <- trace
  OUT$na.action <- attr(mf, "na.action")
  OUT$dataClasses <- attr(attr(mf, "terms"), "dataClasses")[varNamesRHS]
  OUT$model <- mf
  OUT$formula <- formula  
  OUT$parNLS <- unlist(NLS$par) ## coef's from NLS model
  OUT$parONLS <- ONLS$par ## coef's from ONLS model
  OUT$x0 <- X0s ## x_0i
  OUT$y0 <- Y0s ## y_0i
  OUT$fittedONLS  <- fittedONLS ## fitted values from ONLS
  OUT$fittedNLS <- fittedNLS ## fitted values from NLS
  OUT$residONLS <- residONLS ## vertical residuals from ONLS
  OUT$residNLS <- residNLS ## vertical residuals from NLS
  OUT$resid_o  <- resid.o ## orthogonal residuals from ONLS 
  OUT$pred <- PRED
  attr(OUT$pred, "name") <- varNamesRHS
  OUT$resp <- RESP
  attr(OUT$resp, "name") <- varNamesLHS 
  OUT$grad <- GRAD
  OUT$QR <- QR      
  OUT$weights <- wts 
  OUT$control <- control  
  class(OUT) <- c("onls", "nls")
  CHECK <- check_o(OUT, plot = FALSE)
  OUT$ortho <- CHECK$Ortho  
  
  OUT
}