# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

#' At-risk-of-poverty threshold
#' 
#' Estimate the at-risk-of-poverty threshold.  The standard definition is to use
#' 60\% of the weighted median equivalized disposable income.
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
#' @param data an optional \code{data.frame}.
#' @param p a numeric vector of values in \eqn{[0,1]} giving the percentages of 
#' the weighted median to be used for the at-risk-of-poverty threshold.
#' @param na.rm a logical indicating whether missing values should be removed.
#' 
#' @return A numeric vector containing the value(s) of the at-risk-of-poverty
#' threshold is returned.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{arpr}}, \code{\link{incMedian}},
#' \code{\link{weightedMedian}}
#' 
#' @references Working group on Statistics on Income and Living Conditions
#' (2004) Common cross-sectional EU indicators based on EU-SILC; the gender pay
#' gap.  \emph{EU-SILC 131-rev/04}, Eurostat.
#' 
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' arpt("eqIncome", weights = "rb050", data = eusilc)
#' 
#' @export

arpt <- function(inc, weights = NULL, sort = NULL, 
                 years = NULL, data = NULL, p = 0.6, na.rm = FALSE) {
  # check 'p' (other arguments are checked in 'incMedian')
#   if(!is.numeric(p) || length(p) == 0 || p[1] < 0 || p[1] > 1) {
#     stop("'p' must be a numeric value in [0,1]")
#   } else p <- p[1]
  p <- checkP(p)
  byP <- length(p) > 1
  if(byP) {
    if(!is.null(years)) {
      stop("breakdown into years not implemented ",
           "for different threshold levels")
    }
    names(p) <- getPLabels(p)  # ensure that result has correct names
  }
  # compute at-risk-of-poverty threshold
  p * incMedian(inc, weights, sort, years, data, na.rm=na.rm)
}
