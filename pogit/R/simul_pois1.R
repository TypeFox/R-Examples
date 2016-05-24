#' Simulated data set
#' 
#' The data set \code{simul_pois1} contains 300 simulated Poisson counts. 
#' 10 regressors are generated, six of them continuous N(0,1)-variables
#' and four binary with \eqn{p(x_i)=0.5}. The regression effects are set to
#' \code{beta = {2,1,0.6,0,0,1.2,0,0,0.4,-0.2,0.3}}. 
#'  
#' @docType data
#' @usage data(simul_pois1)
#' @format A data frame with 300 rows and the following 12 variables: 
#' \describe{
#'  \item{\code{y}}{number of counts for each covariate pattern}
#'  \item{\code{X.0}}{intercept}
#'  \item{\code{X.1}, \code{X.2}, \code{X.3}, \code{X.4}, \code{X.5}, \code{X.6}, \code{X.7}, \code{X.8}, \code{X.9}, \code{X.10}}{covariates}
#' }
#' 
#' @seealso \code{\link{poissonBvs}} 
#' @name simul_pois1
#' @keywords datasets
NULL
