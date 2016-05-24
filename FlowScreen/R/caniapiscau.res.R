#' @name caniapiscau.res
#' @title Screening results for the Caniapiscau River
#' @description Contains the results from \code{\link{metrics.all}} for the full
#'   Caniapiscau River daily flow series. Data set created as indicated below.  This 
#'   data set is used in the example documentation for the \code{\link{screen.frames}}, 
#'   \code{\link{screen.summary}}, and \code{\link{screen.cpts}} functions in order to 
#'   reduce example run times.
#' @docType data
#' @usage data(caniapiscau)
#' @format Formatted as indicated in the documentation for \code{\link{metrics.all}}
#' @source Original flow series from Environment Canada. 2010. EC Data Explorer V1.2.30. \cr
#'   Water Survey of Canada V1.2.30 https://www.ec.gc.ca/rhc-wsc/
#' @examples
#' # Code used produce this data set:
#' \dontrun{
#' data(caniapiscau)
#' caniapiscau.ts <- create.ts(caniapiscau, hyrstart=3)
#' caniapiscau.ts <- subset(caniapiscau.ts, caniapiscau.ts$hyear > 1962)
#' caniapiscau.res <- metrics.all(caniapiscau.ts)
#' }
#' # example use of example subset flow series
#' data(caniapiscau.res)


"caniapiscau.res"