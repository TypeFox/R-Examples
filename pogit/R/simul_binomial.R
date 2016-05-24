#' Simulated data set
#' 
#' The data set \code{simul_binomial} contains simulated binomial data with 
#' 9 binary covariates. The design matrix is built by computing all 2^9 possible 
#' 0/1 combinations. The regression effects are set to 
#' \code{alpha = {-0.5,0.2,-0.15,0.1,-1.1,0,0,1.2,-0.1,0.3}}. 
#' The number of trials \code{N} are simulated from a Poisson distribution with 
#' parameter \eqn{\exp(\alpha)/(1+\exp(\alpha))}*100.
#' 
#' @docType data
#' @usage data(simul_binomial)
#' @format A data frame with 512 rows and the following 12 variables: 
#' \describe{
#'  \item{\code{y}}{number of successes for each covariate pattern}
#'  \item{\code{N}}{number of trials for each covariate pattern}
#'  \item{\code{X.0}}{intercept}
#'  \item{\code{X.1}, \code{X.2}, \code{X.3}, \code{X.4}, \code{X.5}, \code{X.6}, \code{X.7}, \code{X.8}, \code{X.9}}{binary covariates}
#' }
#' 
#' @seealso \code{\link{logitBvs}}
#' @name simul_binomial
#' @keywords datasets
NULL


