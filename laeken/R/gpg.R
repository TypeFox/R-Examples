# ---------------------------------------
# Author: Matthias Templ
#         Vienna University of Technology
# ---------------------------------------

#' Gender pay (wage) gap.
#' 
#' Estimate the gender pay (wage) gap.
#' 
#' The implementation strictly follows the Eurostat definition (with default
#' method \code{"mean"} and alternative method \code{"median"}).  If weights are
#' provided, the weighted mean or weighted median is estimated.
#' 
#' @param inc either a numeric vector giving the equivalized disposable income,
#' or (if \code{data} is not \code{NULL}) a character string, an integer or a
#' logical vector specifying the corresponding column of \code{data}.
#' @param gender either a factor giving the gender, or (if \code{data} is not
#' \code{NULL}) a character string, an integer or a logical vector specifying
#' the corresponding column of \code{data}.
#' @param method a character string specifying the method to be used.  Possible
#' values are \code{"mean"} for the mean, and \code{"median"} for the median.
#' If weights are provided, the weighted mean or weighted median is estimated.
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
#' @return A list of class \code{"gpg"} (which inherits from the class
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
#' computing the confidence interv al(s) (i.e., the confidence level is \eqn{1 -
#' }\code{alpha}), or \code{NULL}.
#' @returnItem years a numeric vector containing the different years of the
#' survey.
#' @returnItem strata a character vector containing the different domains of the
#' breakdown.
#' 
#' @author Matthias Templ and Alexander Haider, using code for breaking down
#' estimation by Andreas Alfons
#' 
#' @seealso \code{\link{variance}}, \code{\link{qsr}}, \code{\link{gini}}
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
#' data(ses)
#' 
#' # overall value with mean
#' gpg("earningsHour", gender = "sex", weigths = "weights", 
#'     data = ses)
#' 
#' # overall value with median
#' gpg("earningsHour", gender = "sex", weigths = "weights", 
#'     data = ses, method = "median")
#' 
#' # values by education with mean
#' gpg("earningsHour", gender = "sex", weigths = "weights", 
#'     breakdown = "education", data = ses)
#' 
#' # values by education with median
#' gpg("earningsHour", gender = "sex", weigths = "weights", 
#'     breakdown = "education", data = ses, method = "median")
#' 
#' @export

gpg <- function(inc, gender = NULL, method = c("mean", "median"), 
                weights = NULL, sort = NULL, years = NULL, breakdown = NULL, 
                design = NULL, cluster = NULL, data = NULL, var = NULL, 
                alpha = 0.05, na.rm = FALSE, ...) {
  ## initializations
  if(is.null(gender)) stop("'gender' must be supplied")
  byYear <- !is.null(years)
  byStratum <- !is.null(breakdown)
  if(!is.null(data)) {
    inc <- data[, inc]
    gender <- data[, gender]
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
  method <- match.arg(method)
  if(!is.factor(gender)) stop("'gender' must be a factor.")
  if(length(levels(gender)) != 2) stop("'gender' must have exactly two levels")
  if(!all(levels(gender) == c("female", "male"))) {
    gender <- factor(gender, labels=c("female","male"))
    warning("The levels of gender were internally recoded - your first level has to correspond to females")
  }
  if(!is.null(years)) {
    if(!is.factor(years)) stop("'years' should be a factor")
    nage <- length(levels(years))
    if(n > 12) warning(paste("Too small sample sizes may occur by using ", n," age classes"))
  }
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
  # GPG by year (if requested)
  if(byYear) {
    ys <- sort(unique(years))  # unique years
    gp <- function(y, inc, weights, sort, years, na.rm) {
      i <- years == y
      genderGap(inc[i], gender[i], method, weights[i], sort[i], na.rm=na.rm)
    }
    value <- sapply(ys, gp, inc=inc, weights=weights, 
                    sort=sort, years=years, na.rm=na.rm)
    names(value) <- ys  # use years as names
  } else {
    ys <- NULL
    value <- genderGap(inc, gender, method, weights, sort, na.rm=na.rm)
  }
  # GPG by stratum (if requested)
  if(byStratum) {
    gpR <- function(i, inc, weights, sort, na.rm) {
      genderGap(inc[i], gender[i], method, weights[i], sort[i], na.rm=na.rm)
    }
    valueByStratum <- aggregate(1:n, 
                                if(byYear) list(year=years, stratum=breakdown) 
                                else list(stratum=breakdown), 
                                gpR, inc=inc, weights=weights, 
                                sort=sort, na.rm=na.rm)
    names(valueByStratum)[ncol(valueByStratum)] <- "value"
    rs <- levels(breakdown)  # unique strata
  } else valueByStratum <- rs <- NULL
  ## create object of class "gpg"
  res <- constructGpg(value=value, 
                      valueByStratum=valueByStratum, 
                      years=ys, strata=rs)
  # variance estimation (if requested)
  if(!is.null(var)) {
    res <- variance(inc, weights, years, breakdown, design, cluster, 
                    indicator=res, alpha=alpha, na.rm=na.rm, type=var, 
                    gender=gender, method=method, ...)
  }
  ## return result
  return(res)
}

## workhorse
genderGap <- function(x, gend, method = 'mean', weights = NULL, 
                      sort = NULL, na.rm = FALSE) {
  
  if(is.null(gend)) stop("'gender' must be supplied")
  
  # initializations
  if(isTRUE(na.rm)){
    indices <- !is.na(x)
    x <- x[indices]
    gend <- gend[indices]
    if(!is.null(weights)) weights <- weights[indices]
    if(!is.null(sort)) sort <- sort[indices]
  } else if(any(is.na(x))) return(NA)
  
  male <- levels(gend)[1]
  female <- levels(gend)[2]
  
  
  if(is.null(weights)) weights <- rep.int(1, length(x))  # equal weights
  
  incgendmale <- x[gend=="male"]
  incgendmaleWeights <- weights[gend=="male"]
  incgendfemale <- x[gend=="female"]
  incgendfemaleWeights <- weights[gend=="female"]
  
  
  if(method == 'mean') {
    wM <- weighted.mean(x=incgendmale, w=incgendmaleWeights)
    wF <- weighted.mean(x=incgendfemale, w=incgendfemaleWeights)
    return((wM - wF) / wM)
    
  } else {
    wM <- weightedMedian(incgendmale, incgendmaleWeights)
    wF <- weightedMedian(incgendfemale, incgendfemaleWeights)
    return((wM - wF)/wM)
  }
}
