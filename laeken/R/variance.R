# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

#' Variance and confidence intervals of indicators on social exclusion and
#' poverty
#' 
#' Compute variance and confidence interval estimates of indicators on social
#' exclusion and poverty.
#' 
#' This is a wrapper function for computing variance and confidence interval
#' estimates of indicators on social exclusion and poverty.
#' 
#' @param inc either a numeric vector giving the equivalized disposable income,
#' or (if \code{data} is not \code{NULL}) a character string, an integer or a
#' logical vector specifying the corresponding column of \code{data}.
#' @param weights optional; either a numeric vector giving the personal sample
#' weights, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.
#' @param years optional; either a numeric vector giving the different years of
#' the survey, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.  If supplied, values are computed for each year.
#' @param breakdown optional; either a numeric vector giving different domains,
#' or (if \code{data} is not \code{NULL}) a character string, an integer or a
#' logical vector specifying the corresponding column of \code{data}.  If
#' supplied, the values for each domain are computed in addition to the overall
#' value.
#' @param design optional; either an integer vector or factor giving different
#' strata for stratified sampling designs, or (if \code{data} is not
#' \code{NULL}) a character string, an integer or a logical vector specifying
#' the corresponding column of \code{data}.
#' @param cluster optional; either an integer vector or factor giving different
#' clusters for cluster sampling designs, or (if \code{data} is not
#' \code{NULL}) a character string, an integer or a logical vector specifying
#' the corresponding column of \code{data}.
#' @param data an optional \code{data.frame}.
#' @param indicator an object inheriting from the class \code{"indicator"} that
#' contains the point estimates of the indicator (see \code{\link{arpr}},
#' \code{\link{qsr}}, \code{\link{rmpg}} or \code{\link{gini}}).
#' @param alpha a numeric value giving the significance level to be used for
#' computing the confidence interval(s) (i.e., the confidence level is \eqn{1 -
#' }\code{alpha}), or \code{NULL}.
#' @param na.rm a logical indicating whether missing values should be removed.
#' @param type a character string specifying the type of variance estimation to
#' be used.  Currently, only \code{"bootstrap"} is implemented for variance
#' estimation based on bootstrap resampling (see \code{\link{bootVar}}).
#' @param gender either a numeric vector giving the gender, or (if \code{data}
#' is not \code{NULL}) a character string, an integer or a logical vector
#' specifying the corresponding column of \code{data}.
#' @param method a character string specifying the method to be used (only for 
#' \code{\link{gpg}}).  Possible values are \code{"mean"} for the mean, and 
#' \code{"median"} for the median.  If weights are provided, the weighted mean 
#' or weighted median is estimated.
#' @param \dots additional arguments to be passed to \code{\link{bootVar}}.
#' 
#' @return An object of the same class as \code{indicator} is returned.  See
#' \code{\link{arpr}}, \code{\link{qsr}}, \code{\link{rmpg}} or
#' \code{\link{gini}} for details on the components.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{bootVar}}, \code{\link{arpr}}, \code{\link{qsr}},
#' \code{\link{rmpg}}, \code{\link{gini}}
#' 
#' @references 
#' A. Alfons and M. Templ (2013) Estimation of Social Exclusion Indicators 
#' from Complex Surveys: The \R Package \pkg{laeken}.  \emph{Journal of 
#' Statistical Software}, \bold{54}(15), 1--25.  URL 
#' \url{http://www.jstatsoft.org/v54/i15/}
#' 
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' a <- arpr("eqIncome", weights = "rb050", data = eusilc)
#' 
#' ## naive bootstrap
#' variance("eqIncome", weights = "rb050", design = "db040", 
#'     data = eusilc, indicator = a, R = 50, 
#'     bootType = "naive", seed = 123)
#' 
#' ## bootstrap with calibration
#' variance("eqIncome", weights = "rb050", design = "db040", 
#'     data = eusilc, indicator = a, R = 50, 
#'     X = calibVars(eusilc$db040), seed = 123)
#' 
#' @export

variance <- function(inc, weights = NULL, years = NULL, breakdown = NULL, 
                     design = NULL, cluster = NULL, data = NULL, indicator, 
                     alpha = 0.05, na.rm = FALSE, type = "bootstrap", 
                     gender = NULL, method = NULL, ...) {
  # initializations
  type <- match.arg(type)
  # call function corresponding to 'type'
  switch(type,
         bootstrap = bootVar(inc, weights, years, breakdown, design, cluster, 
                             data, indicator, alpha=alpha, na.rm=na.rm, 
                             gender=gender, method=method, ...))
}

