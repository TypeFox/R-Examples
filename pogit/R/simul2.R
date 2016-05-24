#' Simulated data set
#' 
#' The simulated data set \code{simul2} considers a situation with clustered 
#' observations and four binary covariates in both sub-models of the Pogit
#' model, i.e. 
#' \code{X} = \code{W}. The respective design matrix is
#' built by computing all 2^4 possible 0/1 combinations and one observation 
#' is generated for each covariate pattern. C=50 clusters are built containing 
#' one unit with each of the resulting 16 covariate patterns, i.e. a total of 
#' I=800 units. The regression effects are set to \code{beta = {0.75,0.1,0.1,0,0}}
#' in the Poisson and to \code{alpha = {2.2,-0.3,0,-0.3,0}} in the logit model. 
#' Random intercepts in both sub-models are simulated from a normal distribution
#' with standard deviations \eqn{\theta_\beta}=\code{0.1} and 
#' \eqn{\theta_\alpha}=\code{0.3}. Additionally to the main study sample, 
#' validation data are available for each covariate pattern and cluster. 
#' For details concerning the simulation setup, see Dvorzak and Wagner (2016). 
#' 
#' @docType data
#' @usage data(simul2)
#' @format A data frame with 800 rows and the following 10 variables: 
#' \describe{
#'  \item{\code{y}}{number of observed counts for each covariate pattern 
#'    in each cluster}
#'  \item{\code{E}}{total exposure times for each unit}
#'  \item{\code{cID}}{cluster ID for each unit}
#'  \item{\code{X.0}}{intercept}
#'  \item{\code{X.1}, \code{X.2}, \code{X.3}, \code{X.4}}{binary covariates}
#'  \item{\code{v}}{number of reported cases for each covariate pattern in each
#'    cluster in the validation sample}
#'  \item{\code{m}}{number of true cases subject to the fallible reporting
#'    process (sample size of validation data)}
#' }
#' 
#' @source Dvorzak, M. and Wagner, H. (2016). Sparse Bayesian modelling
#'  of underreported count data. \emph{Statistical Modelling}, \strong{16}(1),
#'  24 - 46, \url{http://dx.doi.org/10.1177/1471082x15588398}. 
#' 
#' @seealso \code{\link{pogitBvs}}
#' @name simul2
#' @keywords datasets
NULL