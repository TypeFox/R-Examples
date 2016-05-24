# aliases meanWt meanWt.default meanWt.dataObj varWt varWt.default
# varWt.dataObj covWt covWt.default covWt.matrix covWt.data.frame
# covWt.dataObj corWt corWt.default corWt.matrix corWt.data.frame
# corWt.dataObj

#' @name weighted_estimators
#' @rdname weighted_estimators
#' @title Weighted mean, variance, covariance matrix and correlation matrix
#' @description Compute mean, variance, covariance matrix and correlation matrix, taking
#' into account sample weights.
#' \itemize{
#' \item \code{meanWt}: a simple wrapper that calls \code{mean(x, na.rm=na.rm)} if
#' \code{weights} is missing and \code{weighted.mean(x, w=weights,
#' na.rm=na.rm)} otherwise. Implemented methods for this generic are:
#' \itemize{
#' \item \code{meanWt.default(x, weights, na.rm=TRUE, ...)}
#' \item \code{meanWt.dataObj(x, vars, na.rm=TRUE, ...)}
#' }
#' \item \code{varWt}: calls \code{var(x, na.rm=na.rm)} if \code{weights} is missing.
#' Implemented methods for this generic are:
#' \itemize{
#' \item \code{varWt.default(x, weights, na.rm=TRUE, ...)}
#' \item \code{varWt.dataObj(x, vars, na.rm=TRUE, ...)}
#' }
#' \item \code{covWt} and \code{covWt}: always remove missing values pairwise and call
#' \code{cov} and \code{cor}, respectively, if \code{weights} is missing.
#' Implemented methods for these generics are:
#' \itemize{
#' \item \code{covWt.default(x, y, weights, ...)}
#' \item \code{covWt.matrix(x, weights, ...)}
#' \item \code{covWt.data.frame(x, weights, ...) }
#' \item \code{covWt.dataObj(x, vars, ...)}
#' \item \code{corWt.default(x, y, weights, ...)}
#' \item \code{corWt.matrix(x, weights, ...)}
#' \item \code{corWt.data.frame(x, weights, ...)}
#' \item \code{corWt.dataObj(x, vars, ...)}
#' }
#' }
#' The additional parameters are now described:
#' \itemize{
#' \item y: a numeric vector.  If missing, this defaults to \code{x}.
#' \item vars: a character vector of variable names that should be used for the
#' calculation.
#' \item na.rm: a logical indicating whether any \code{NA} or \code{NaN} values
#' should be removed from \code{x} before computation.  Note that the default
#' is \code{TRUE}.
#' \item weights: an optional numeric vector containing sample weights.
#' }
#'
#' @param x for \code{meanWt} and \code{varWt}, a numeric vector or an object
#' of class \code{\linkS4class{dataObj}}. For \code{covWt} and \code{corWt}, a
#' numeric vector, matrix, \code{data.frame} or \code{\linkS4class{dataObj}}.
#' In case of a \code{\linkS4class{dataObj}}, weights are automatically used
#' from the S4-object itself.
#' @param \dots for the generic functions \code{covWt} and \code{corWt},
#' additional arguments to be passed to methods.  Additional arguments not
#' included in the definition of the methods are ignored.
#' @return For \code{meanWt}, the (weighted) mean.
#'
#' For \code{varWt}, the (weighted) variance.
#'
#' For \code{covWt}, the (weighted) covariance matrix or, for the default
#' method, the (weighted) covariance.
#'
#' For \code{corWt}, the (weighted) correlation matrix or, for the default
#' method, the (weighted) correlation coefficient.
#' @author Stefan Kraft and Andreas Alfons
#' @seealso \code{\link{mean}}, \code{\link[stats]{weighted.mean}},
#' \code{\link[stats:cor]{var}}, \code{\link[stats:cor]{cov}},
#' \code{\link[stats:cor]{cor}}
#' @keywords univar multivariate array
#' @note \code{meanWt}, \code{varWt}, \code{covWt} and \code{corWt} all make use of
#' slot \code{weights} of the input object if the \code{dataObj}-method is
#' used.
NULL

#' @rdname weighted_estimators
#' @examples
#' data(eusilcS)
#' meanWt(eusilcS$netIncome, weights=eusilcS$rb050)
#' sqrt(varWt(eusilcS$netIncome, weights=eusilcS$rb050))
#'
#' # dataObj-methods
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' meanWt(inp, vars="netIncome")
#' sqrt(varWt(inp, vars="netIncome"))
#' corWt(inp, vars=c("age", "netIncome"))
#' covWt(inp, vars=c("age", "netIncome"))
#' @export
meanWt <- function(x, ...) UseMethod("meanWt")

#' @export
meanWt.default <- function(x, weights, na.rm=TRUE, ...) {
  na.rm <- isTRUE(na.rm)
  if(missing(weights)) mean(x, na.rm=na.rm)
  else weighted.mean(x, w=weights, na.rm=na.rm)
}

#' @export
meanWt.dataObj <- function(x, vars, na.rm=TRUE, ...) {
  dat <- x@data
  if ( is.null(dat) ) {
    return(NULL)
  } else {
    if ( length(vars) > 1 ) {
      stop("only one variable can be specified!\n")
    }
    ii <- match(vars, colnames(dat))
    if ( any(is.na(ii)) ) {
      stop("please provide valid variables that exist in the input object!\n")
    }
    tmpdat <- dat[[vars]]
    if ( !is.null(x@weight) ) {
      return(meanWt.default(tmpdat, weights=dat[[x@weight]], na.rm=na.rm))
    } else {
      return(meanWt.default(tmpdat, na.rm=na.rm))
    }
  }
}

## weighted variance
#' @rdname weighted_estimators
#' @export
varWt <- function(x, ...) UseMethod("varWt")

#' @export
varWt.default <- function(x, weights, na.rm=TRUE, ...) {
  na.rm <- isTRUE(na.rm)
  if(missing(weights)) var(x, na.rm=na.rm)
  else {
    x <- as.numeric(x)
    weights <- as.numeric(weights)
    if(length(weights) != length(x)) {
      stop("'weights' must have the same length as 'x'")
    }
    if(na.rm) {
      select <- !is.na(x)
      x <- x[select]
      weights <- weights[select]
    }
    if(length(x) <= 1 || sum(weights > 0) <= 1) NA
    else sum((x - meanWt(x, weights))^2 * weights) / (sum(weights) - 1)
  }
}

#' @export
varWt.dataObj <- function(x, vars, na.rm=TRUE, ...) {
  dat <- x@data
  if ( is.null(dat) ) {
    return(NULL)
  } else {
    if ( length(vars) > 1 ) {
      stop("only one variable can be specified!\n")
    }
    ii <- match(vars, colnames(dat))
    if ( any(is.na(ii)) ) {
      stop("please provide valid variables that exist in the input object!\n")
    }
    tmpdat <- dat[[vars]]
    if ( !is.null(x@weight) ) {
      return(varWt.default(tmpdat, weights=dat[[x@weight]], na.rm=na.rm))
    } else {
      return(varWt.default(tmpdat, na.rm=na.rm))
    }
  }
}

## weighted covariance matrix
## generic function
#' @rdname weighted_estimators
#' @export
covWt <- function(x, ...) UseMethod("covWt")


## default method
#' @export
covWt.default <- function(x, y, weights, ...) {
  if(missing(y)) y <- x
  else if(length(x) != length(y)) {
    stop("'x' and 'y' must have the same length")
  }
  if(missing(weights)) cov(x, y, use = "complete.obs")
  else {
    if(length(weights) != length(x)) {
      stop("'weights' must have the same length as 'x' and 'y'")
    }
    select  <- !is.na(x) & !is.na(y)
    x <- x[select]
    y <- y[select]
    weights <- weights[select]
    sum((x-meanWt(x, weights)) * (y-meanWt(y, weights)) * weights) /
      (sum(weights)-1)
  }
}


### method for matrices
#' @export
covWt.matrix <- function(x, weights, ...) {
  if(missing(weights)) cov(x, use = "pairwise.complete.obs")
  else {
    if(length(weights) != nrow(x)) {
      stop("length of 'weights' must equal the number of rows in 'x'")
    }
    center <- apply(x, 2, meanWt, weights=weights)
    x <- sweep(x, 2, center, check.margin = FALSE)
    crossprodWt(x, weights)
  }
}


### method for data.frames
#' @export
covWt.data.frame <- function(x, weights, ...) covWt(as.matrix(x), weights)


#' @export
covWt.dataObj <- function(x, vars, ...) {
  dat <- x@data
  if ( is.null(dat) ) {
    return(NULL)
  } else {
    ii <- match(vars, colnames(dat))
    if ( any(is.na(ii)) ) {
      stop("please provide valid variables that exist in the input object!\n")
    }
    tmpdat <- dat[,vars,with=F]
    if ( !is.null(x@weight) ) {
      return(covWt.matrix(as.matrix(tmpdat), weights=dat[[x@weight]]))
    } else {
      return(covWt.matrix(as.matrix(tmpdat)))
    }
  }
}


## weighted correlation matrix

### generic function
#' @rdname weighted_estimators
#' @export
corWt <- function(x, ...) UseMethod("corWt")

### default method
#' @export
corWt.default <- function(x, y, weights, ...) {
  if(missing(y)) y <- x
  else if(length(x) != length(y)) {
    stop("'x' and 'y' must have the same length")
  }
  if(missing(weights)) cor(x, y, use = "complete.obs")
  else covWt(x, y, weights) / sqrt(varWt(x, weights) * varWt(y, weights))
}


#### method for matrices
#' @export
corWt.matrix <- function(x, weights, ...) {
  if(missing(weights)) cor(x, use = "pairwise.complete.obs")
  else {
    if(length(weights) != nrow(x)) {
      stop("length of 'weights' must equal the number of rows in 'x'")
    }
    cen <- apply(x, 2, meanWt, weights=weights)
    sc <- sqrt(apply(x, 2, varWt, weights=weights))
    x <- scale(x, center=cen, scale=sc)
    crossprodWt(x, weights)
  }
}


### method for data.frames
#' @export
corWt.data.frame <- function(x, weights, ...) corWt(as.matrix(x), weights)


### method for objects of class "dataObj"
#' @export
corWt.dataObj <- function(x, vars, ...) {
  dat <- x@data
  if ( is.null(dat) ) {
    return(NULL)
  } else {
    ii <- match(vars, colnames(dat))
    if ( any(is.na(ii)) ) {
      stop("please provide valid variables that exist in the input object!\n")
    }
    tmpdat <- dat[,vars,with=F]
    if ( !is.null(x@weight) ) {
      return(corWt.matrix(as.matrix(tmpdat), weights=dat[[x@weight]]))
    } else {
      return(corWt.matrix(as.matrix(tmpdat)))
    }
  }
}


### weighted cross product
### designed for internal use, hence no error handling
### TODO: more efficient solution
crossprodWt <- function(x, weights) {
  ci <- 1:ncol(x)
  sapply(ci, function(j) sapply(ci, function(i) {
    select <- !is.na(x[, i]) & !is.na(x[, j])
    xi <- x[select, i]
    xj <- x[select, j]
    w <- weights[select]
    sum(xi*xj*w) / (sum(w)-1)
  }))
}
