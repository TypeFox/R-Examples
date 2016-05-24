#' The SPU and aSPU tests for multiple traits - single SNP association in generalized estimating equations.
#'
#' It gives p-values of the GEESPU tests and GEEaSPU test.
#'
#' @param traits trait matrix. The row for individuals and the column for traits.
#'
#' @param geno A matrix of genetic information.
#'
#' @param Z covariates.
#'
#' @param model Use "gaussian" for a quantitative trait, and use "binomial" for a binary trait.
#'
#' @param gamma power used in GEEaSPU test. A vector of the powers.
#'
#' @param n.sim number of simulations.
#'
#' @param corstr a character string specifying the correlation structure. The following are permitted: "independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", "AR-M" and "unstructured" 
#'
#' @return p-values for the GEE-SPU and GEE-aSPU test.
#'
#' @author Junghi Kim, Wei Pan and Il-Youp Kwak
#'
#' @references
#'
#' Yiwei Zhang, Zhiyuan Xu, Xiaotong Shen, Wei Pan (2014)
#' Testing for association with multiple traits in generalized estimation equations, with application to neuroimaging data.
#' Neuroimage. 96:309-25 
#'
#' @examples
#'
#' traits <- matrix(rnorm(100*5, 0,1), ncol=5)
#' Z <- rnorm(100, 2, 0.5)
#' geno <- rbinom(100, 2, 0.5)
#' out <- GEEaSPU(traits, geno, Z = NULL, model = "gaussian", 
#'		  gamma = c(1:8,Inf), n.sim = 100)
#'

GEEaSPU <- function(traits, geno, 
	  Z = NULL, model = c("binomial", "gaussian"), 
	  gamma = c(1:8,Inf), n.sim=1000, corstr = "independence") {

   Score <- GEE(traits = traits, geno = geno, Z = Z,family = model, corstr = corstr)
   U <- Score$U
   V <- Score$Cov
   out <- GEEspu.score(U = U, V = V, gamma = gamma, B = n.sim)
   out
}

