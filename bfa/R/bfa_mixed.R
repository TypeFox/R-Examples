#' @include main.R
NULL
#' Initialize and fit a mixed-scale Gaussian factor model, with probit specifications for the discrete
#' margins
#'
#' This function performs a specified number of MCMC iterations and
#' returns an object containing summary statistics from the MCMC samples
#' as well as the actual samples if keep.scores or keep.loadings are \code{TRUE}.
#' Default behavior is to save only the loadings. 
#' 
#' Additional parameters:
#' \itemize{
#' \item loadings.var: Factor loading prior variance
#' \item tau.a, tau.b: Gamma hyperparameters (scale=1/b) for factor precisions (if factor.scales=T). Default is a=b=1 (MV t w/ df=2)
#' \item rho.a, rho.b: Beta hyperparameters for point mass prior
#' \item sigma2.a, sigma2.b: Gamma hyperparameters for error precisions (for numeric variables)
#' \item gdp.alpha, gdp.beta: GDP prior parameters
#' }
#' 
#' @param x A formula or bfa object. 
#' @param data The data if x is a formula
#' @param num.factor Number of factors
#' @param restrict A matrix or list giving restrictions on factor loadings. A matrix should be the 
#' same size as the loadings matrix. Acceptable values are 0 (identically 0), 1 (unrestricted), 
#' or 2 (strictly positive). List elements should be character vectors of the form c('variable',1, ">0")
#' where 'variable' is the manifest variable, 1 is the factor, and ">0" is the restriction. Acceptable
#' restrictions are ">0" or "0".
#' @param normal.dist A character vector specifying which variables should be treated as observed 
#' Gaussian. Defaults to all numeric variables if x is a formula.
#' @param center.data Center data (for continuous variables - recommended!)
#' @param scale.data  Scale data (for continuous variables)
#' @param nsim Number of iterations past burn-in
#' @param nburn Number of initial (burn-in) iterations to discard
#' @param thin Keep every thin'th MCMC sample (i.e. save nsim/thin samples)
#' @param print.status How often to print status messages to console
#' @param keep.scores Save samples of factor scores
#' @param loading.prior Specify GDP ("gdp", default) point mass ("pointmass") or normal priors ("normal") 
#' @param factor.scales Include a separate scale parameter for each factor
#' @param px Use parameter expansion for ordinal variables (recommended)
#' @param coda Create \code{mcmc} objects to allow use of functions from the 
#' \code{coda} package: "all" for loadings and scores, "loadings" or "scores" for one or the
#' other, or "none" for neither
#' @param coda.scale Put the loadings on the correlation scale when creating \code{mcmc} objects
#' @param imh.iter Iterations used to build IMH proposal
#' @param imh.burn Burn-in before collecting samples used to build IMH proposal (total burn-in is nburn+imh.iter+imh.burn)
#' @param ... Prior parameters and other (experimental) arguments (see details)
#' @return A \code{bfa} object with posterior samples.
#' @export

bfa_mixed <- function(x, data=NULL, num.factor=1, restrict=NA, normal.dist=NA, 
                      center.data=TRUE, scale.data=TRUE, nsim=10000, nburn=1000, thin=1,
                      print.status=500, keep.scores=FALSE,
                      loading.prior=c("gdp", "pointmass", "normal"), 
                      factor.scales=FALSE, px=TRUE,
                      coda="loadings", coda.scale=TRUE, imh.iter=500,
                      imh.burn=500, ...) {
  keep.loadings = TRUE
  loading.prior = match.arg(loading.prior)
  m = .bfa(x, data=data, num.factor=num.factor, restrict=restrict, normal.dist=normal.dist, 
           center.data=center.data, scale.data=scale.data, nsim=nsim, nburn=nburn, thin=thin,
           print.status=print.status, keep.scores=keep.scores, keep.loadings=keep.loadings,
           loading.prior=loading.prior, factor.scales=factor.scales, px=px, coda=coda, 
           coda.scale=coda.scale, imh=TRUE, imh.iter=imh.iter, imh.burn=imh.burn, ...)
  attr(m, "type") <- "mixed"
  if(nsim==10000 & nburn==1000 & thin==1) warning("Default MCMC settings used; assess convergence! (see ?get_coda)")
  return(m)
}
