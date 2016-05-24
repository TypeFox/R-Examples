#' Radon data
#' 
#' Radon measurements of 919 owner-occupied homes in 85 counties of Minnesota.
#' 
#' @usage data(radon)
#' @docType data
#' @name radon
#' @keywords datasets
#' @format A data frame with 919 observations on the following 5 variables:
#' \describe{
#' \item{log.radon}{Radon measurement (in log pCi/L, i.e., log picoCurie per liter)}
#' \item{basement}{Indicator for the level of the home at which the radon measurement
#'   was taken - 0 = basement, 1 = first floor.}
#' \item{uranium}{Average county-level soil uranium content.}
#' \item{county}{County ID.}
#' \item{county.name}{County name - a factor.}
#' }
#' @source 
#' \url{http://www.stat.columbia.edu/~gelman/arm/software/}
#' @references
#'   Price, P. N., Nero, A. V. and Gelman, A. (1996) Bayesian prediction of mean 
#'   indoor radon concentrations for Minnesota counties. \emph{Health Physics}.
#'   \bold{71}(6), 922--936.
#'   
#'   Gelman, A. and Hill, J. (2007) \emph{Data analysis using regression and 
#'   multilevel/hierarchical models}. Cambridge University Press.
NULL