nlsModel <- function (form, data, start, wts, upper = NULL) 
{
    thisEnv <- environment()
    env <- new.env(hash = TRUE, parent = environment(form))
    for (i in names(data)) assign(i, data[[i]], envir = env)
    ind <- as.list(start)
    parLength <- 0
    for (i in names(ind)) {
        temp <- start[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp, envir = env)
        ind[[i]] <- parLength + seq_along(start[[i]])
        parLength <- parLength + length(start[[i]])
    }
    getPars.noVarying <- function() unlist(setNames(lapply(names(ind), 
        get, envir = env), names(ind)))
    getPars <- getPars.noVarying
    internalPars <- getPars()
    if (!is.null(upper)) upper <- rep_len(upper, parLength)
    useParams <- rep(TRUE, parLength)
    lhs <- eval(form[[2L]], envir = env)
    rhs <- eval(form[[3L]], envir = env)
    .swts <- if (!missing(wts) && length(wts)) sqrt(wts) else rep_len(1, length(rhs))
    assign(".swts", .swts, envir = env)
    resid <- .swts * (lhs - rhs)
    dev <- sum(resid^2)
    if (is.null(attr(rhs, "gradient"))) {
        getRHS.noVarying <- function() {
            if (is.null(upper)) 
                numericDeriv(form[[3L]], names(ind), env)
            else numericDeriv(form[[3L]], names(ind), env, ifelse(internalPars < 
                upper, 1, -1))
        }
        getRHS <- getRHS.noVarying
        rhs <- getRHS()
    }
    else {
        getRHS.noVarying <- function() eval(form[[3L]], envir = env)
        getRHS <- getRHS.noVarying
    }
    dimGrad <- dim(attr(rhs, "gradient"))
    marg <- length(dimGrad)
    if (marg > 0L) {
        gradSetArgs <- vector("list", marg + 1L)
        for (i in 2L:marg) gradSetArgs[[i]] <- rep(TRUE, dimGrad[i - 
            1])
        useParams <- rep(TRUE, dimGrad[marg])
    }
    else {
        gradSetArgs <- vector("list", 2L)
        useParams <- rep(TRUE, length(attr(rhs, "gradient")))
    }
    npar <- length(useParams)
    gradSetArgs[[1L]] <- (~attr(ans, "gradient"))[[2L]]
    gradCall <- switch(length(gradSetArgs) - 1L, call("[", gradSetArgs[[1L]], 
        gradSetArgs[[2L]], drop = FALSE), call("[", gradSetArgs[[1L]], 
        gradSetArgs[[2L]], gradSetArgs[[2L]], drop = FALSE), 
        call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]], 
            gradSetArgs[[3L]], drop = FALSE), call("[", gradSetArgs[[1L]], 
            gradSetArgs[[2L]], gradSetArgs[[2L]], gradSetArgs[[3L]], 
            gradSetArgs[[4L]], drop = FALSE))
    getRHS.varying <- function() {
        ans <- getRHS.noVarying()
        attr(ans, "gradient") <- eval(gradCall)
        ans
    }
    QR <- qr(.swts * attr(rhs, "gradient"))
    qrDim <- min(dim(QR$qr))
    if (QR$rank < qrDim) stop("singular gradient matrix at initial parameter estimates")
    
    getPars.varying <- function() unlist(setNames(lapply(names(ind), 
        get, envir = env), names(ind)))[useParams]
    setPars.noVarying <- function(newPars) {
        assign("internalPars", newPars, envir = thisEnv)
        for (i in names(ind)) assign(i, unname(newPars[ind[[i]]]), 
            envir = env)
    }
    setPars.varying <- function(newPars) {
        internalPars[useParams] <- newPars
        for (i in names(ind)) assign(i, unname(internalPars[ind[[i]]]), 
            envir = env)
    }
    setPars <- setPars.noVarying
    on.exit(remove(i, data, parLength, start, temp, m))
    m <- list(resid = function() resid, fitted = function() rhs, 
        formula = function() form, deviance = function() dev, 
        lhs = function() lhs, gradient = function() .swts * attr(rhs, 
            "gradient"), conv = function() {
            if (npar == 0) return(0)
            rr <- qr.qty(QR, resid)
            sqrt(sum(rr[1L:npar]^2)/sum(rr[-(1L:npar)]^2))
        }, incr = function() qr.coef(QR, resid), setVarying = function(vary = rep(TRUE, 
            length(useParams))) {
            assign("useParams", if (is.character(vary)) {
                temp <- logical(length(useParams))
                temp[unlist(ind[vary])] <- TRUE
                temp
            } else if (is.logical(vary) && length(vary) != length(useParams)) stop("setVarying : 'vary' length must match length of parameters") else {
                vary
            }, envir = thisEnv)
            gradCall[[length(gradCall) - 1L]] <<- useParams
            if (all(useParams)) {
                assign("setPars", setPars.noVarying, envir = thisEnv)
                assign("getPars", getPars.noVarying, envir = thisEnv)
                assign("getRHS", getRHS.noVarying, envir = thisEnv)
                assign("npar", length(useParams), envir = thisEnv)
            } else {
                assign("setPars", setPars.varying, envir = thisEnv)
                assign("getPars", getPars.varying, envir = thisEnv)
                assign("getRHS", getRHS.varying, envir = thisEnv)
                assign("npar", length(seq_along(useParams)[useParams]), 
                  envir = thisEnv)
            }
        }, setPars = function(newPars) {
            setPars(newPars)
            assign("resid", .swts * (lhs - assign("rhs", getRHS(), 
                envir = thisEnv)), envir = thisEnv)
            assign("dev", sum(resid^2), envir = thisEnv)
            assign("QR", qr(.swts * attr(rhs, "gradient")), envir = thisEnv)
            return(QR$rank < min(dim(QR$qr)))
        }, getPars = function() getPars(), getAllPars = function() getPars(), 
        getEnv = function() env, trace = function() {
            cat(format(dev), ": ", format(getPars()))
            cat("\n")
        }, Rmat = function() qr.R(QR), predict = function(newdata = list(), 
            qr = FALSE) eval(form[[3L]], as.list(newdata), env))
    class(m) <- "nlsModel"
    m
}

nlsLM <- function (formula, data = parent.frame(), start, jac = NULL, 
                  algorithm = "LM", control = nls.lm.control(), 
                  lower = NULL, upper = NULL, trace = FALSE, subset, weights, na.action, 
                  model = FALSE, ...) 
{
  formula <- as.formula(formula)
  if (!is.list(data) && !is.environment(data)) stop("'data' must be a list or an environment")
  mf <- match.call()
  varNames <- all.vars(formula)
  
  if (length(formula) == 2L) {
    formula[[3L]] <- formula[[2L]]
    formula[[2L]] <- 0
  }

  form2 <- formula
  form2[[2L]] <- 0
  varNamesRHS <- all.vars(form2)
  mWeights <- missing(weights)
  
  ## if trace = TRUE, set nls.lm.control$nprint = 1
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
  } else names(start)

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
      else stop("parameters without starting value in 'data': ", paste(nnn, collapse = ", "))
    }
  } else {
    if (length(pnames) && any((np <- sapply(pnames, lenVar)) == -1)) {
      message("fitting parameters ", paste(sQuote(pnames[np == -1]), collapse = ", "), " without any variables")
      n <- integer()
    } else stop("no parameters to fit")
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
      wts <- if (mWeights) rep(1, n) else eval(substitute(weights), data, environment(formula))
    } else {
      mf$formula <- as.formula(paste("~", paste(varNames[varIndex], collapse = "+")), env = environment(formula))
      mf$start <- mf$control <- mf$algorithm <- mf$trace <- mf$model <- NULL
      mf$lower <- mf$upper <- NULL
      mf[[1L]] <- as.name("model.frame")
      mf <- eval.parent(mf)
      n <- nrow(mf)
      mf <- as.list(mf)
      wts <- if (!mWeights) model.weights(mf) else rep(1, n)
    }
    if (any(wts < 0 | is.na(wts))) stop("missing or negative weights not allowed")
  } else {
    varIndex <- logical()
    mf <- list(0)
    wts <- numeric()
  }
  
  if (missing(start)) start <- getInitial(formula, mf)
  for (var in varNames[!varIndex]) mf[[var]] <- eval(as.name(var), data, env)
  varNamesRHS <- varNamesRHS[varNamesRHS %in% varNames[varIndex]]

  ## added nls.lm from 'minpack.lm' package
  mf <- c(mf, start)
  lhs <- eval(formula[[2L]], envir = mf)
  m <- match(names(start), names(mf)) 
  .swts <- if (!missing(wts) && length(wts)) sqrt(wts)
  
  FCT <- function(par) {    
    mf[m] <- par
    rhs <- eval(formula[[3L]], envir = mf, environment(formula))
    res <- lhs - rhs
    res <- .swts * res    
    res
  }
  
  NLS <- nls.lm(par = start, fn = FCT, jac = jac, control = control, lower = lower, upper = upper, ...)
  ## previous versions before boundaries were included
  # NLS <- nls.lm(par = start, fn = FCT, jac = jac, control = control, ...)
  start <- NLS$par
        
  ##pass optimized parameters to 'nlsModel'
  m <- nlsModel(formula, mf, start, wts)
    
  ## => internal 'nls' iterations by 'C_nls_iter' switched off
  ## use 'nls.lm' output for convergence info
  if (NLS$info %in% c(1, 2, 3, 4)) isConv <- TRUE else isConv <- FALSE
  finIter <- NLS$niter   
  finTol <- nls.lm.control()$ftol
  convInfo <- list(isConv = isConv, finIter = finIter, finTol = finTol, 
                    stopCode = NLS$info, stopMessage = NLS$message)
  nls.out <- list(m = m, convInfo = convInfo, data = substitute(data), call = match.call())

  nls.out$call$algorithm <- algorithm
  ## need to use '$tol' parameter from nls.control to make 'predict.nls' work
  nls.out$call$control <- nls.control()
  nls.out$call$trace <- FALSE
  nls.out$na.action <- attr(mf, "na.action")
  nls.out$dataClasses <- attr(attr(mf, "terms"), "dataClasses")[varNamesRHS]
  if (model) 
    nls.out$model <- mf
  if (!mWeights) 
    nls.out$weights <- wts
  nls.out$control <- control
  class(nls.out) <- "nls"
  nls.out
}
