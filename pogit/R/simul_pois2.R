#' Simulated data set
#' 
#' The same simulation setup is used as in \code{\link{simul_pois1}} but considers 
#' clustered observations. 10 regressors are generated, six of them continuous 
#' N(0,1)-variables and four binary with \eqn{p(x_i)=0.5}. 
#' The regression effects are set to \code{beta = {2,1,0.6,0,0,1.2,0,0,0.4,-0.2,0.3}}.
#' To simulate clustering, it is assumed that each of 
#' C=10 clusters is formed of 30 subjects and 10 random intercepts are generated 
#' from a normal distribution with zero mean and standard deviation 
#' \eqn{\theta} = 0.1. 
#'  
#' @docType data
#' @usage data(simul_pois2)
#' @format A data frame with 300 rows and the following 12 variables: 
#' \describe{
#'  \item{\code{y}}{number of counts for each covariate pattern in each cluster}
#'  \item{\code{cID}}{cluster ID of each count}
#'  \item{\code{X.0}}{intercept}
#'  \item{\code{X.1}, \code{X.2}, \code{X.3}, \code{X.4}, \code{X.5}, \code{X.6}, \code{X.7}, \code{X.8}, \code{X.9}, \code{X.10}}{covariates}
#' }
#' 
#' @seealso \code{\link{simul_pois1}}, \code{\link{poissonBvs}}
#' @name simul_pois2
#' @keywords datasets
NULL

