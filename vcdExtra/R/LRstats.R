# fixed buglet when deviance() returns a null
# fixed bug: residual df calculated incorrectly
#  but this now depends on objects having a df.residual component
#  TRUE for lm, glm, polr, negbin objects

# made generic, adding a glmlist method

LRstats <- function(object, ...) {
	UseMethod("LRstats")
}

LRstats.glmlist <- function(object, ..., saturated = NULL, sortby=NULL)
{
    ns <- sapply(object, function(x) length(x$residuals))
    if (any(ns != ns[1L])) 
        stop("models were not all fitted to the same size of dataset")
    nmodels <- length(object)
    if (nmodels == 1) 
        return(LRstats.default(object[[1L]], saturated=saturated))
    
    rval <- lapply(object, LRstats.default, saturated=saturated)
    rval <- do.call(rbind, rval)    
		if (!is.null(sortby)) {
			rval <- rval[order(rval[,sortby], decreasing=TRUE),]
			}
		rval
}

# could just do LRstats.loglmlist <- LRstats.glmlist
LRstats.loglmlist <- function(object, ..., saturated = NULL, sortby=NULL)
{
	ns <- sapply(object, function(x) length(x$residuals))
	if (any(ns != ns[1L])) 
		stop("models were not all fitted to the same size of dataset")
	nmodels <- length(object)
	if (nmodels == 1) 
		return(LRstats.default(object[[1L]], saturated=saturated))
	
	rval <- lapply(object, LRstats.default, saturated=saturated)
	rval <- do.call(rbind, rval)    
	if (!is.null(sortby)) {
		rval <- rval[order(rval[,sortby], decreasing=TRUE),]
	}
	rval
}

LRstats.default <- function(object, ..., saturated = NULL, sortby=NULL)
{
  ## interface methods for logLik() and nobs()
  ## - use S4 methods if loaded
  ## - use residuals() if nobs() is not available
  logLik0 <- if("stats4" %in% loadedNamespaces()) stats4::logLik else logLik
  nobs0   <- function(x, ...) {
    nobs1 <- if("stats4" %in% loadedNamespaces()) stats4::nobs else nobs
    nobs2 <- function(x, ...) NROW(residuals(x, ...))
    rval <- try(nobs1(x, ...), silent = TRUE)
    if(inherits(rval, "try-error") | is.null(rval)) rval <- nobs2(x, ...)
    return(rval)
  }
  dof <- function(x) {
  	if (inherits(x, "loglm")) {
  		rval <- x$df 
  		} else {
  		rval <- try(x$df.residual, silent=TRUE)
  		}
  	if (inherits(rval, "try-error") || is.null(rval)) stop(paste("Can't determine residual df for a", class(x), "object"))
  	rval
  	}

  ## collect all objects
  objects <- list(object, ...)
  nmodels <- length(objects)
  
  ## check sample sizes
  ns <- sapply(objects, nobs0)
  if(any(ns != ns[1L])) stop("models were not all fitted to the same size of dataset")

  ## extract log-likelihood and df (number of parameters)
  ll <- lapply(objects, logLik0)
  par <- as.numeric(sapply(ll, function(x) attr(x, "df")))
	df <- as.numeric(sapply(objects, function(x) dof(x)))
  ll <- sapply(ll, as.numeric)
  
  ## compute saturated reference value (use 0 if deviance is not available)
  if(is.null(saturated)) {
    dev <- try(sapply(objects, deviance), silent = TRUE)
    if(inherits(dev, "try-error") || any(sapply(dev, is.null))) {
      saturated <- 0
    } else {
      saturated <- ll + dev/2
    }
  }

  ## setup ANOVA-style matrix
  rval <- matrix(rep(NA, 5 * nmodels), ncol = 5)
  colnames(rval) <- c("AIC", "BIC", "LR Chisq", "Df", "Pr(>Chisq)")
  rownames(rval) <- as.character(sapply(match.call(), deparse)[-1L])[1:nmodels]
  rval[,1] <- -2 * ll + 2 * par
  rval[,2] <- -2 * ll + log(ns) * par
  rval[,3] <- -2 * (ll - saturated)
  rval[,4] <- df
  rval[,5] <- pchisq(rval[,3], df, lower.tail = FALSE)

	if (!is.null(sortby)) {
		rval <- rval[order(rval[,sortby], decreasing=TRUE),]
	}

  ## return
  structure(as.data.frame(rval), heading = "Likelihood summary table:", class = c("anova", "data.frame"))
}
