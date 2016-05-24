#' R package for fast and exact (restricted) likelihood ratio tests for mixed and additive models.
#'
#' \code{RLRsim} implements fast simulation-based exact tests for variance components in mixed and additive models for 
#' conditionally Gaussian responses -- i.e., tests for questions like: 
#' \itemize{
#' \item is the variance of my random intercept significantly different from 0?
#' \item is this smooth effect significantly nonlinear?
#' \item is this smooth effect significantly different from a constant effect?}
#' The convenience functions \code{\link{exactRLRT}} and \code{\link{exactLRT}}
#'  can deal with fitted models from packages \pkg{lme4, nlme, gamm4, SemiPar} and 
#'  from \pkg{mgcv}'s \code{gamm()}-function.
#' Workhorse functions  \code{\link{LRTSim}} and  \code{\link{RLRTSim}} 
#' accept design matrices as inputs directly and can thus be used more generally
#' to generate exact critical values for the corresponding 
#' (restricted) likelihood ratio tests.\cr\cr
#' The theory behind these tests was first developed in:\cr
#' Crainiceanu, C. and Ruppert, D. (2004) 
#' \href{http://people.orie.cornell.edu/~davidr/papers/asymptoticpaper2.pdf}{Likelihood ratio tests in
#' linear mixed models with one variance component}, \emph{Journal of the Royal
#' Statistical Society: Series B}, \bold{66}, 165--185.\cr\cr
#' Power analyses and sensitivity studies for \pkg{RLRsim} can be found in:\cr
#' Scheipl, F., Greven, S. and Kuechenhoff, H. (2008) 
#' \href{http://dx.doi.org/10.1016/j.csda.2007.10.022}{Size and power of tests
#' for a zero random effect variance or polynomial regression in additive and
#' linear mixed models}.  \emph{Computational Statistics and Data Analysis},
#' \bold{52}(7), 3283--3299.
#' 
#' 
#' 
#' @name RLRsim-package
#' @aliases RLRsim-package RLRsim
#' @docType package
#' @author Fabian Scheipl (\email{fabian.scheipl@@stat.uni-muenchen.de}), 
#'   Ben Bolker
#' @keywords package
NULL



