waldtest <- function(object, ...) {
  UseMethod("waldtest")
}

waldtest.formula <- function(object, ..., data = list()) {
  object <- if(length(data) < 1) eval(call("lm", formula = as.formula(deparse(substitute(object))),
    environment(object)))
  else eval(call("lm", formula = as.formula(deparse(substitute(object))),
    data = as.name(deparse(substitute(data))), environment(data)))
  waldtest.lm(object, ...)
}

waldtest.default <- function(object, ..., vcov = NULL, test = c("Chisq", "F"), name = NULL)
{
  ## methods needed:
  ## - terms()
  ## - update()
  ## - formula()
  ## - nobs() or residuals() -> only for determining number of observations
  ## - df.residual() or logLik
  ## - coef() -> needs to be named, matching names in terms() and vcov()
  ## - vcov(), potentially user-supplied

  ## use S4 methods if loaded
  coef0   <- if("stats4" %in% loadedNamespaces()) stats4::coef   else coef
  logLik0 <- if("stats4" %in% loadedNamespaces()) stats4::logLik else logLik
  update0 <- if("stats4" %in% loadedNamespaces()) stats4::update else update
  nobs0   <- function(x, ...) {
    nobs1 <- if("stats4" %in% loadedNamespaces()) stats4::nobs else nobs
    nobs2 <- function(x, ...) NROW(residuals(x, ...))
    rval <- try(nobs1(x, ...), silent = TRUE)
    if(inherits(rval, "try-error") | is.null(rval)) rval <- nobs2(x, ...)
    return(rval)
  }
  vcov0   <- if(!is.null(vcov)) vcov else {
    if("stats4" %in% loadedNamespaces()) stats4::vcov else stats::vcov
  }
  df.residual0 <- function(x) {
    df <- try(df.residual(x), silent = TRUE)
    if(inherits(df, "try-error") | is.null(df)) df <- try(nobs0(x) - attr(logLik0(x), "df"), silent = TRUE)
    if(inherits(df, "try-error") | is.null(df)) df <- try(nobs0(x) - length(as.vector(coef0(x))), silent = TRUE)
    if(inherits(df, "try-error")) df <- NULL
    return(df)
  }

  ## model class
  cls <- class(object)[1]

  ## convenience functions:
  ## 1. extracts term labels
  tlab <- function(x) {
    tt <- try(terms(x), silent = TRUE)
    if(inherits(tt, "try-error")) "" else attr(tt, "term.labels")
  }
  ## 2. extracts model name
  if(is.null(name)) name <- function(x) {
    rval <- try(formula(x), silent = TRUE)
    if(inherits(rval, "try-error") | is.null(rval)) rval <- try(x$call, silent = TRUE)
    if(inherits(rval, "try-error") | is.null(rval)) return(NULL) else return(paste(deparse(rval), collapse="\n"))
  }
  ## 3. compute an updated model object
  modelUpdate <- function(fm, update) {
    ## if `update' is numeric or character, then assume that the 
    ## corresponding variables (i.e., terms) are redundant (i.e., should be omitted)
    if(is.numeric(update)) {
      ## sanity checking of numeric update specification
      if(any(update < 1)) {
        warning("for numeric model specifications all values have to be >=1")
	update <- abs(update)[abs(update) > 0]
      }
      if(any(update > length(tlab(fm)))) {
        warning(paste("more terms specified than existent in the model:",
	        paste(as.character(update[update > length(tlab(fm))]), collapse = ", ")))
	update <- update[update <= length(tlab(fm))]
      }
      ## finally turn numeric into character update specification
      update <- tlab(fm)[update]
    }
    if(is.character(update)) {
      ## sanity checking of character update specification
      if(!all(update %in% tlab(fm))) {
        warning(paste("terms specified that are not in the model:",
	        paste(dQuote(update[!(update %in% tlab(fm))]), collapse = ", ")))
        update <- update[update %in% tlab(fm)]
      }
      if(length(update) < 1) stop("empty model specification")  
      ## finally turn character into formula update specification       
      update <- as.formula(paste(". ~ . -", paste(update, collapse = " - ")))
    }
    if(inherits(update, "formula")) update <- update0(fm, update)
    if(!inherits(update, cls)) stop(paste("original model was of class \"", cls,
      "\", updated model is of class \"", class(update)[1], "\"", sep = ""))
    return(update)
  }
  ## 4. compare two fitted model objects
  modelCompare <- function(fm, fm.up, vfun = NULL) {
    q <- length(coef0(fm)) - length(coef0(fm.up))

    if(q > 0) {
      fm0 <- fm.up
      fm1 <- fm
    } else {
      fm0 <- fm
      fm1 <- fm.up
    }
    k <- length(coef0(fm1))
    n <- nobs0(fm1)

    ## determine omitted variables
    if(!all(tlab(fm0) %in% tlab(fm1))) stop("models are not nested")
    ovar <- which(!(names(coef0(fm1)) %in% names(coef0(fm0))))

    ## get covariance matrix estimate
    vc <- if(is.null(vfun)) vcov(fm1)
          else if(is.function(vfun)) vfun(fm1)
	  else vfun

    ## compute Chisq statistic
    stat <- t(coef0(fm1)[ovar]) %*% solve(vc[ovar,ovar]) %*% coef0(fm1)[ovar]
    return(c(-q, stat))
  }

  ## recursively fit all objects (if necessary)
  objects <- list(object, ...)
  nmodels <- length(objects)
  if(nmodels < 2) {
    objects <- c(objects, . ~ 1)
    nmodels <- 2
  }
  
  # remember which models are already fitted and which are described
  # by an update mechanism
  no.update <- sapply(objects, function(obj) inherits(obj, cls))
  
  ## updating
  for(i in 2:nmodels) objects[[i]] <- modelUpdate(objects[[i-1]], objects[[i]])

  ## check responses
  getresponse <- function(x) {
    tt <- try(terms(x), silent = TRUE)
    if(inherits(tt, "try-error")) "" else deparse(tt[[2]])
  }
  responses <- as.character(lapply(objects, getresponse))
  sameresp <- responses == responses[1]
  if(!all(sameresp)) {
    objects <- objects[sameresp]
    warning("models with response ", deparse(responses[!sameresp]),
	    " removed because response differs from ", "model 1")
  }

  ## check sample sizes
  ns <- sapply(objects, nobs0)
  if(any(ns != ns[1])) {
    for(i in 2:nmodels) {
      if(ns[1] != ns[i]) {
        if(no.update[i]) stop("models were not all fitted to the same size of dataset")
	  else {
	    commonobs <- row.names(model.frame(objects[[i]])) %in% row.names(model.frame(objects[[i-1]]))
	    objects[[i]] <- eval(substitute(update(objects[[i]], subset = commonobs),
	      list(commonobs = commonobs)))
	    if(nobs0(objects[[i]]) != ns[1]) stop("models could not be fitted to the same size of dataset")
	  }
      }
    }
  }

  ## check vcov0
  if(nmodels > 2 && !is.null(vcov0) && !is.function(vcov0))
    stop("to compare more than 2 models `vcov' needs to be a function")

  ## setup ANOVA matrix
  test <- match.arg(test)
  rval <- matrix(rep(NA, 4 * nmodels), ncol = 4)
  colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
  rownames(rval) <- 1:nmodels
  rval[,1] <- as.numeric(sapply(objects, df.residual0))
  for(i in 2:nmodels) rval[i, 2:3] <- modelCompare(objects[[i-1]], objects[[i]], vfun = vcov0)
  if(test == "Chisq") {
    rval[,4] <- pchisq(rval[,3], round(abs(rval[,2])), lower.tail = FALSE)
  } else {
    df <- rval[,1]
    for(i in 2:nmodels) if(rval[i,2] < 0) df[i] <- rval[i-1,1]
    rval[,3] <- rval[,3]/abs(rval[,2])
    rval[,4] <- pf(rval[,3], abs(rval[,2]), df, lower.tail = FALSE)
  }

  variables <- lapply(objects, name)
  if(any(sapply(variables, is.null))) variables <- lapply(match.call()[-1L], deparse)[1L:nmodels]
  title <- "Wald test\n"
  topnote <- paste("Model ", format(1:nmodels),": ", variables, sep="", collapse="\n")

  structure(as.data.frame(rval), heading = c(title, topnote),
	    class = c("anova", "data.frame"))
}

waldtest.lm <- function(object, ..., test = c("F", "Chisq"))
{
  waldtest.default(object, ..., test = match.arg(test))
}
