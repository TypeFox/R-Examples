# ---------------------------------------
# Author: Matthias Templ
#         Vienna University of Technology
# ---------------------------------------

#' Proportion of an alternative distribution
#'
#' Estimate the proportion of an alternative distribution.
#'
#' If weights are provided, the weighted proportion is estimated.
#'
#' @param bin either a factor vector giving the values,
#' or (if \code{data} is not \code{NULL}) a character string, an integer or a
#' logical vector specifying the corresponding column of \code{data}.
#' @param weights optional; either a numeric vector giving the personal sample
#' weights, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.
#' @param sort optional; either a numeric vector giving the personal IDs to be
#' used as tie-breakers for sorting, or (if \code{data} is not \code{NULL}) a
#' character string, an integer or a logical vector specifying the corresponding
#' column of \code{data}.
#' @param years optional; either a numeric vector giving the different years of
#' the survey, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.  If supplied, values are computed for each year.
#' @param breakdown optional; either a numeric vector giving different domains,
#' or (if \code{data} is not \code{NULL}) a character string, an integer or a
#' logical vector specifying the corresponding column of \code{data}.  If
#' supplied, the values for each domain are computed in addition to the overall
#' value.
#' @param design optional and only used if \code{var} is not \code{NULL}; either
#' an integer vector or factor giving different domains for stratified sampling
#' designs, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.
#' @param cluster optional and only used if \code{var} is not \code{NULL};
#' either an integer vector or factor giving different clusters for cluster
#' sampling designs, or (if \code{data} is not \code{NULL}) a character string,
#' an integer or a logical vector specifying the corresponding column of
#' \code{data}.
#' @param data an optional \code{data.frame}.
#' @param var a character string specifying the type of variance estimation to
#' be used, or \code{NULL} to omit variance estimation.  See
#' \code{\link{variance}} for possible values.
#' @param alpha numeric; if \code{var} is not \code{NULL}, this gives the
#' significance level to be used for computing the confidence interval (i.e.,
#' the confidence level is \eqn{1 - }\code{alpha}).
#' @param na.rm a logical indicating whether missing values should be removed.
#' @param \dots if \code{var} is not \code{NULL}, additional arguments to be
#' passed to \code{\link{variance}}.
#'
#' @return A list of class \code{"prop"} (which inherits from the class
#' \code{"indicator"}) with the following components:
#' @returnItem value a numeric vector containing the overall value(s).
#' @returnItem valueByStratum a \code{data.frame} containing the values by
#' domain, or \code{NULL}.
#' @returnItem varMethod a character string specifying the type of variance
#' estimation used, or \code{NULL} if variance estimation was omitted.
#' @returnItem var a numeric vector containing the variance estimate(s), or
#' \code{NULL}.
#' @returnItem varByStratum a \code{data.frame} containing the variance
#' estimates by domain, or \code{NULL}.
#' @returnItem ci a numeric vector or matrix containing the lower and upper
#' endpoints of the confidence interval(s), or \code{NULL}.
#' @returnItem ciByStratum a \code{data.frame} containing the lower and upper
#' endpoints of the confidence intervals by domain, or \code{NULL}.
#' @returnItem alpha a numeric value giving the significance level used for
#' computing the confidence interval(s) (i.e., the confidence level is \eqn{1 -
#' }\code{alpha}), or \code{NULL}.
#' @returnItem years a numeric vector containing the different years of the
#' survey.
#' @returnItem strata a character vector containing the different domains of the
#' breakdown.
#'
#' @author Matthias Templ, using code for breaking down
#' estimation by Andreas Alfons
#'
#' @seealso \code{\link{variance}}
#'
#' @references
#' A. Alfons and M. Templ (2013) Estimation of Social Exclusion Indicators
#' from Complex Surveys: The \R Package \pkg{laeken}.  \emph{Journal of
#' Statistical Software}, \bold{54}(15), 1--25.  URL
#' \url{http://www.jstatsoft.org/v54/i15/}
#'
#' Working group on Statistics on Income and Living Conditions (2004)
#' Common cross-sectional EU indicators based on EU-SILC; the gender
#' pay gap.  \emph{EU-SILC 131-rev/04}, Eurostat, Luxembourg.
#'
#' @keywords survey
#'
#' @examples
#' data(eusilc)
#'
#' # overall value
#' prop("rb090", weights = "rb050", data = eusilc)
#'
#' # values by region
#' prop("rb090", weights = "rb050",
#'     breakdown = "db040", data = eusilc)
#'
#' @export

prop <- function(bin, weights = NULL, sort = NULL, years = NULL,
                 breakdown = NULL, design = NULL, cluster = NULL,
                 data = NULL, var = NULL, alpha = 0.05,
                 na.rm = FALSE, ...) {
  ## initializations
  byYear <- !is.null(years)
  byStratum <- !is.null(breakdown)
  if(!is.null(data)) {
    bin <- data[, bin]
    if(!is.null(weights)) weights <- data[, weights]
    if(!is.null(sort)) sort <- data[, sort]
    if(byYear) years <- data[, years]
    if(byStratum) breakdown <- data[, breakdown]
    if(!is.null(var)) {
      if(!is.null(design)) design <- data[, design]
      if(!is.null(cluster)) cluster <- data[, cluster]
    }
  }
  # check vectors
  if(!is.factor(bin)) stop("'bin' must be a vector of class factor")
  if(length(levels(bin)) != 2) stop(paste("'bin' has not exactly 2 levels"))
  n <- length(bin)
  if(is.null(weights)) weights <- weights <- rep.int(1, n)
  else if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
  if(!is.null(sort) && !is.vector(sort) && !is.ordered(sort)) {
    stop("'sort' must be a vector or ordered factor")
  }
  if(byYear && !is.numeric(years)) {
    stop("'years' must be a numeric vector")
  }
  if(byStratum) {
    if(!is.vector(breakdown) && !is.factor(breakdown)) {
      stop("'breakdown' must be a vector or factor")
    } else breakdown <- as.factor(breakdown)
  }
  if(is.null(data)) {  # check vector lengths
    if(length(weights) != n) {
      stop("'weights' must have the same length as 'x'")
    }
    if(!is.null(sort) && length(sort) != n) {
      stop("'sort' must have the same length as 'x'")
    }
    if(byYear && length(years) != n) {
      stop("'years' must have the same length as 'x'")
    }
    if(byStratum && length(breakdown) != n) {
      stop("'breakdown' must have the same length as 'x'")
    }
  }
  ## computations
  # prop by year (if requested)
  if(byYear) {
    ys <- sort(unique(years))  # unique years
    gc <- function(y, bin, weights, sort, years, na.rm) {
      i <- years == y
      propCoeff(bin[i], weights[i], sort[i], na.rm=na.rm)
    }
    value <- sapply(ys, gc, bin=bin, weights=weights,
                    sort=sort, years=years, na.rm=na.rm)
    names(value) <- ys  # use years as names
  } else {
    ys <- NULL
    value <- propCoeff(bin, weights, sort, na.rm=na.rm)
  }
  # prop by stratum (if requested)
  if(byStratum) {
    gcR <- function(i, bin, weights, sort, na.rm) {
      propCoeff(bin[i], weights[i], sort[i], na.rm=na.rm)
    }
    valueByStratum <- aggregate(1:n,
                                if(byYear) list(year=years, stratum=breakdown)
                                else list(stratum=breakdown),
                                gcR, bin=bin, weights=weights,
                                sort=sort, na.rm=na.rm)
    names(valueByStratum)[ncol(valueByStratum)] <- "value"
    rs <- levels(breakdown)  # unique strata
  } else valueByStratum <- rs <- NULL
  ## create object of class "qsr"
  res <- constructProp(value=value,
                       valueByStratum=valueByStratum,
                       years=ys, strata=rs)
  # variance estimation (if requested)
  if(!is.null(var)) {
    bin <- as.numeric(as.integer(bin) - 1)
    res <- variance(bin, weights, years, breakdown, design, cluster,
                    indicator=res, alpha=alpha, na.rm=na.rm, type=var, ...)
  }
  ## return result
  return(res)
}

## workhorse
propCoeff <- function(x, weights = NULL, sort = NULL, na.rm = FALSE) {
  # initializations
  if(isTRUE(na.rm)){
    indices <- !is.na(x)
    x <- x[indices]
    if(!is.null(weights)) weights <- weights[indices]
    if(!is.null(sort)) sort <- sort[indices]
  } else if(any(is.na(x))) return(NA)
  # sort values and weights
#  order <- if(is.null(sort)) order(x) else order(x, sort)
#  x <- x[order]  # order values
  if(is.null(weights)) weights <- rep.int(1, length(x))  # equal weights
#  else weights <- weights[order]  # order weights
  ## calculations
  ## bin to numeric
  x <- as.integer(x)
  weightedMean(x-1, weights)
}
