# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

#' Relative median at-risk-of-poverty gap
#' 
#' Estimate the relative median at-risk-of-poverty gap, which is defined as the
#' relative difference between the median equivalized disposable income of
#' persons below the at-risk-of-poverty threshold and the at-risk-of-poverty
#' threshold itself (expressed as a percentage of the at-risk-of-poverty
#' threshold).
#' 
#' The implementation strictly follows the Eurostat definition.
#' 
#' @param inc either a numeric vector giving the equivalized disposable income,
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
#' value.  Note that the same (overall) threshold is used for all domains.
#' @param design optional and only used if \code{var} is not \code{NULL}; either
#' an integer vector or factor giving different strata for stratified sampling
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
#' @return A list of class \code{"rmpg"} (which inherits from the class
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
#' @returnItem threshold a numeric vector containing the at-risk-of-poverty
#' threshold(s).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{arpt}}, \code{\link{variance}}
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
#' rmpg("eqIncome", weights = "rb050", data = eusilc)
#' 
#' # values by region
#' rmpg("eqIncome", weights = "rb050", 
#'     breakdown = "db040", data = eusilc)
#' 
#' @export

rmpg <- function(inc, weights = NULL, sort = NULL, years = NULL, 
                 breakdown = NULL, design = NULL, cluster = NULL, 
                 data = NULL, var = NULL, alpha = 0.05, 
                 na.rm = FALSE, ...) {
  ## initializations
  byYear <- !is.null(years)
  byStratum <- !is.null(breakdown)
  if(!is.null(data)) {
    inc <- data[, inc]
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
  if(!is.numeric(inc)) stop("'inc' must be a numeric vector")
  n <- length(inc)
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
  if(byYear) {  # RMPG by year
    ys <- sort(unique(years))
    ts <- arpt(inc, weights, sort, years, na.rm=na.rm)  # thresholds
    rg <- function(y, t, inc, weights, sort, years, na.rm) {
      i <- years == y
      relativeGap(inc[i], weights[i], sort[i], t, na.rm=na.rm)
    }
    value <- mapply(rg, y=ys, t=ts, 
                    MoreArgs=list(inc=inc, weights=weights, sort=sort, 
                                  years=years, na.rm=na.rm))
    names(value) <- ys  # use years as names
    if(byStratum) {
      rg1 <- function(i, inc, weights, sort, years, ts, na.rm) {
        y <- years[i[1]]
        t <- ts[as.character(y)]
        relativeGap(inc[i], weights[i], sort[i], t, na.rm=na.rm)
      }
      valueByStratum <- aggregate(1:n, list(year=years, stratum=breakdown), 
                                  rg1, inc=inc, weights=weights, sort=sort, 
                                  years=years, ts=ts, na.rm=na.rm)
      names(valueByStratum)[3] <- "value"
    } else valueByStratum <- NULL
  } else {  # RMPG for only one year
    ys <- NULL
    ts <- arpt(inc, weights, sort, na.rm=na.rm)  # threshold
    value <- relativeGap(inc, weights, sort, ts, na.rm=na.rm)
    if(byStratum) {
      rg2 <- function(i, inc, weights, sort, ts, na.rm) {
        relativeGap(inc[i], weights[i], sort[i], ts, na.rm=na.rm)
      }
      valueByStratum <- aggregate(1:n, list(stratum=breakdown), 
                                  rg2, inc=inc, weights=weights, sort=sort, 
                                  ts=ts, na.rm=na.rm)
      names(valueByStratum)[2] <- "value"
    } else valueByStratum <- NULL
  }
  rs <- levels(breakdown)  # unique strata (also works if 'breakdown' is NULL)
  ## create object of class "arpr"
  res <- constructRmpg(value=value, valueByStratum=valueByStratum, 
                       years=ys, strata=rs, threshold=ts)
  # variance estimation (if requested)
  if(!is.null(var)) {
    res <- variance(inc, weights, years, breakdown, design, cluster, 
                    indicator=res, alpha=alpha, na.rm=na.rm, type=var, ...)
  }
  ## return result
  return(res)
}

## workhorse
relativeGap <- function(x, weights = NULL, sort = NULL, 
                        threshold, na.rm = FALSE) {
  ## initializations
  if(is.null(weights)) weights <- rep.int(1, length(x))  # equal weights
  if(isTRUE(na.rm)){
    indices <- !is.na(x)
    x <- x[indices]
    if(!is.null(weights)) weights <- weights[indices]
    if(!is.null(sort)) sort <- sort[indices]
  } else if(any(is.na(x))) return(NA)
  if(length(x) == 0) return(NA)
  # preparations
  isPoor <- x < threshold  # individuals below threshold
  x <- x[isPoor]
  if(!is.null(weights)) weights <- weights[isPoor]
  if(!is.null(sort)) sort <- sort[isPoor]
  # calculations
  medianPoor <- incMedian(x, weights, sort)
  (threshold - medianPoor) * 100 / threshold
}
