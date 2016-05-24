lrtest <- function(object, ...) {
  UseMethod("lrtest")
}

lrtest.formula <- function(object, ..., data = list()) {
  object <- if(length(data) < 1) eval(call("lm", formula = as.formula(deparse(substitute(object))),
    environment(object)))
  else eval(call("lm", formula = as.formula(deparse(substitute(object))),
    data = as.name(deparse(substitute(data))), environment(data)))
  lrtest.default(object, ...)
}

lrtest.default <- function(object, ..., name = NULL)
{
  ## methods needed:
  ## - logLik()

  ## - terms()  -> for updating only
  ## - update()

  ## - formula() -> for determining `name' (can be user-supplied)
  ## - nobs() or residuals() -> for determining number of observations

  ## use S4 methods if loaded
  logLik0 <- if("stats4" %in% loadedNamespaces()) stats4::logLik else logLik
  update0 <- if("stats4" %in% loadedNamespaces()) stats4::update else update
  nobs0   <- function(x, ...) {
    nobs1 <- if("stats4" %in% loadedNamespaces()) stats4::nobs else nobs
    nobs2 <- function(x, ...) NROW(residuals(x, ...))
    rval <- try(nobs1(x, ...), silent = TRUE)
    if(inherits(rval, "try-error") | is.null(rval)) rval <- nobs2(x, ...)
    return(rval)
  }

  ## model class
  cls <- class(object)[1]

  ## convenience functions:
  ## 1. extracts term labels
  tlab <- function(x) attr(terms(x), "term.labels")
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
    if(!inherits(update, cls)) warning(paste("original model was of class \"", cls,
      "\", updated model is of class \"", class(update)[1], "\"", sep = ""))
    return(update)
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

  ## setup ANOVA matrix
  rval <- matrix(rep(NA, 5 * nmodels), ncol = 5)
  colnames(rval) <- c("#Df", "LogLik", "Df", "Chisq", "Pr(>Chisq)")
  rownames(rval) <- 1:nmodels
  
  logL <- lapply(objects, logLik0)
  rval[,1] <- as.numeric(sapply(logL, function(x) attr(x, "df")))  
  rval[,2] <- sapply(logL, as.numeric)
  rval[2:nmodels, 3] <- rval[2:nmodels, 1] - rval[1:(nmodels-1), 1]
  rval[2:nmodels, 4] <- 2 * abs(rval[2:nmodels, 2] - rval[1:(nmodels-1), 2])
  rval[,5] <- pchisq(rval[,4], round(abs(rval[,3])), lower.tail = FALSE)

  variables <- lapply(objects, name)
  if(any(sapply(variables, is.null))) variables <- lapply(match.call()[-1L], deparse)[1L:nmodels]
  title <- "Likelihood ratio test\n"
  topnote <- paste("Model ", format(1:nmodels),": ", variables, sep="", collapse="\n")

  structure(as.data.frame(rval), heading = c(title, topnote),
	    class = c("anova", "data.frame"))
}
