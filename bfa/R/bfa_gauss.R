#' @include main.R
NULL
#' Initialize and fit a Gaussian factor model
#'
#' Performs MCMC model fitting for a Gaussian factor model.
#' 
#' Additional parameters:
#' \itemize{
#' \item loadings.var: Factor loading prior variance
#' \item tau.a, tau.b: Gamma hyperparameters (scale=1/b) for factor precisions (if factor.scales=T). Default is a=b=1 (MV t w/ df=2)
#' \item rho.a, rho.b: Beta hyperparameters for point mass prior
#' \item sigma2.a, sigma2.b: Gamma hyperparameters for error precisions
#' \item gdp.alpha, gdp.beta: GDP prior parameters
#' }
#' 
#' @param x A formula or bfa object. 
#' @param data The data (if x is a formula)
#' @param num.factor Number of factors
#' @param restrict A matrix or list giving identifiability restrictions on factor loadings. 
#' A matrix should be the same size as the loadings matrix. Acceptable values are 0 (identically 0),
#' 1 (unrestricted), 
#' or 2 (strictly positive). List elements should be character vectors of the form c("variable",1, ">0")
#' where 'variable' is the manifest variable, 1 is the factor, and ">0" is the restriction. Acceptable
#' restrictions are ">0" or "0".
#' @param center.data Center data (recommended; this model doesn't include a mean vector!)
#' @param scale.data  Scale data (also recommended)
#' @param nsim Number of iterations past burn-in
#' @param nburn Number of initial (burn-in) iterations to discard
#' @param thin Keep every thin'th MCMC sample (i.e. save nsim/thin samples)
#' @param print.status How often to print status messages to console
#' @param keep.scores Save samples of factor scores
#' @param loading.prior Specify the prior on factor loadings -  generalized double Pareto ("gdp", default),
#' point mass mixtures (mixture of point mass at zero + mean zero normal) ("pointmass") or normal/Gaussian ("normal") 
#' @param factor.scales Include a shared precision parameter for each column of the factor 
#' loadings matrix. See details for setting hyperprior parameters. This is implemented as in PX-FA of
#' Ghosh and Dunson (2009)
#' @param coda Create \code{mcmc} objects to allow use of functions from the 
#' \code{coda} package: "all" for loadings and scores, "loadings" or "scores" for one or the
#' other, or "none" for neither
#' @param ... Prior parameters and other (experimental) arguments (see details)
#' @return A \code{bfa} object with posterior samples.
#' @export

bfa_gauss <- function(x, data=NULL, num.factor=1, restrict=NA, 
                      center.data=TRUE, scale.data=TRUE, nsim=10000, nburn=1000, thin=1,
                      print.status=500, keep.scores=FALSE,
                      loading.prior=c("gdp", "pointmass", "normal"), 
                      factor.scales = TRUE,
                      coda="loadings", ...) {
  keep.loadings=TRUE
  loading.prior = match.arg(loading.prior)
  fr = model.frame(x, data=data, na.action=NULL)
  d = dim(fr)
  normal.dist = rep(1, d[2])
  m = .bfa(x, data=data, num.factor=num.factor, restrict=restrict, normal.dist=normal.dist, 
           center.data=center.data, scale.data=scale.data, nsim=nsim, nburn=nburn, thin=thin,
           print.status=print.status, keep.scores=keep.scores, keep.loadings=keep.loadings,
           loading.prior=loading.prior, factor.scales=factor.scales, px=FALSE, coda=coda, 
           imh=FALSE, ...) 
  attr(m, "type") <- "gauss"
  if(nsim==10000 & nburn==1000 & thin==1) warning("Default MCMC settings used; assess convergence! (see ?get_coda)")
  return(m)         
}
