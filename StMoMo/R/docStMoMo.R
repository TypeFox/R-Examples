#' StMoMo.
#'
#' @name StMoMo
#' @docType package
#' @importFrom gnm gnm
#' @importFrom gnm Mult
#' @importFrom forecast forecast
#' @importFrom stats simulate
#' @importFrom rootSolve multiroot
#' @importFrom fanplot fan
#' @importFrom fields image.plot
#' @importFrom graphics abline frame image.default lines par plot
#' @importFrom grDevices colorRampPalette grey.colors rgb
#' @importFrom stats aggregate as.formula binomial coef cov fitted lm na.omit 
#' @importFrom stats nobs poisson predict qnorm quantile rbinom residuals
#' @importFrom stats rpois runif start
#' @importFrom utils head tail
#' @importFrom RColorBrewer brewer.pal
NULL


# #' England and Wales male mortality data
# #'
# #' Age-specific deaths and exposures for England and Wales from 
# #' the Human Mortality Database.
# #'
# #' \code{EWMaleData} contains deaths and exposures for England and 
# #' Wales males for the period 1961-2011 and for ages 0-100. 
# #' Data taken from the Human Mortality Database on 5 November 2014.
# #'
# #' @format A list  with the following components:
# #' \describe{
# #'   \item{Dxt}{ matrix of deaths data.}
# #'   \item{Ext}{ matrix of exposures data (mid year population estimates).}
# #'   \item{ages}{ vector of ages.}
# #'   \item{years}{ vector of years.}
# #' }
# #' @source Human Mortality Database \url{http://www.mortality.org/}.
# #' @references Human Mortality Database (2014). University of California, 
# #' Berkeley (USA), and Max Planck Institute for Demographic Research (Germany). 
# #' Available at \url{www.mortality.org}.
# "EWMaleData"