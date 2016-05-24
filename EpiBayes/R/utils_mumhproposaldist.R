#' @title 
#' Proposal Distribution for MH Step for \code{mu}
#' 
#' @description 
#' This function caculates realizations from the proposal distribution for \code{mu} in 
#'     the Metropolis-Hastings (MH) step in the hybrid Gibbs sampler used in the 
#'     \code{\link{EpiBayes_ns}} and \code{\link{EpiBayes_s}} functions. Recall,
#'     in the context of disease freedom, \code{mu} was the average disease prevalence
#'     in a given diseased cluster of subjects.
#'
#' @param mumh The proposed or old value of \code{mu}.
#' @param mumh.alpha The old first shape parameter.
#' @param mumh.beta The old second shape parameter.
#' @param mumh.alphastar The new first shape parameter.
#' @param mumh.betastar The new second shape parameter.
#' 
#' @return
#' A realization from the (beta) proposal distribution.
#'
#' @references
#' Norets, A. and Tang, X. \emph{MCMC estimation of a finite beta mixture}. Technical 
#'     report, Princeton University, 2010.
#'
#' @seealso
#' Used in the functions \code{\link{EpiBayes_ns}} and \code{\link{EpiBayes_s}}.
#'
#' @export 
#' @name utils_mumhproposaldist
#' @import stats
utils_mumhproposaldist = function(mumh, mumh.alpha, mumh.beta, mumh.alphastar, mumh.betastar){
						return(dbeta(mumh, mumh.alphastar + mumh.alpha - 1, mumh.betastar + mumh.beta - 1))
					}