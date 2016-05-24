#' @title 
#' Target Distribution for MH Step for \code{mu}
#' 
#' @description 
#' This function caculates realizations from the target distribution for \code{mu} in 
#'     the Metropolis-Hastings (MH) step in the hybrid Gibbs sampler used in the 
#'     \code{\link{EpiBayes_ns}} and \code{\link{EpiBayes_s}} functions. Recall,
#'     in the context of disease freedom, \code{mu} was the average disease prevalence
#'     in a given diseased cluster of subjects.
#'
#' @param mumh The proposed or old value of \code{mu}.
#' @param psi The current MCMC iteration's value for \code{psi} (related to the variability
#'     of the disease prevalences among infected clusters of subjects).
#' @param pi The current MCMC iteration's value for \code{pi} (the subject-level 
#'     prevalences in each cluster).
#' @param z.pi The current MCMC iteration's value for \code{z.pi} (the indicator variables
#'     corresponding to the disease status of each cluster).
#' @param muparm The old first and second shape parameters.
#' 
#' @return
#' A realization from the target distribution.
#'
#' @references
#' Norets, A. and Tang, X. \emph{MCMC estimation of a finite beta mixture}. Technical 
#'     report, Princeton University, 2010.
#'
#' @seealso
#' Used in the functions \code{\link{EpiBayes_ns}} and \code{\link{EpiBayes_s}}.
#' 
#' @export 
#' @name utils_mumhtargetdist
#' @import stats
utils_mumhtargetdist = function(mumh, psi, pi, z.pi, muparm){
						return(prod((1 / (gamma(mumh * psi) * (gamma((1 - mumh) * psi))) * pi^mumh * (1 - pi)^(1 - mumh))^z.pi) * mumh^(muparm[1] - 1) * (1 - mumh)^(muparm[2] - 1))
					}