#'Permutation-based variance component test statistic for longitudinal RNA-seq data
#'
#'This function computes an approximation of the Variance Component test for a
#'mixture of \eqn{\chi^{2}}s using permutations. This is preferable to the
#'asymptotic approximation for small sample sizes
#'
#'
#'@param y a numeric matrix of dim \code{g x n} containing the raw RNAseq counts for g
#'genes from \code{n} samples.
#'
#'@param x a numeric design matrix of dim \code{n x p} containing the \code{p} covariates
#' to be adjusted for
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}.
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the \code{K} variables
#'to be tested
#'
#'@param w a vector of length \code{n} containing the weights for the \code{n}
#'samples.
#'
#'@param Sigma_xi a matrix of size \code{K x K} containing the covariance matrix
#'of the \code{K} random effects.
#'
#'@param n_perm the number of perturbations
#'
#'@return A list with the following elements:\itemize{
#'   \item \code{scores_perm}: a vector of length \code{nperm} containing the
#'   aproximated score test statistics observed for each permutation
#'   \item \code{score_obs}: approximation of the observed score
#'   \item \code{pval}: associated p-value
#' }
#'
#'@seealso \code{\link[CompQuadForm]{davies}}
#'
#'@examples
#'#rm(list=ls())
#'set.seed(123)
#'
#'##generate some fake data
#'########################
#'n <- 100
#'r <- 12
#'t <- matrix(rep(1:3), 4, ncol=1, nrow=r)
#'sigma <- 0.5
#'b0 <- 1
#'
#'#under the null:
#'b1 <- 0
#'#under the alternative:
#'b1 <- 0.7
#'y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
#'y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#'       matrix(rep(y.tilde, n), ncol=n, nrow=r))
#'x <- matrix(1, ncol=1, nrow=r)
#'
#'#run test
#'permTestRes <- vc_test_perm(y, x, phi=t, w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                            indiv=rep(1:4, each=3), n_perm=100)
#'permTestRes$pval
#'
#'@importFrom CompQuadForm davies
#'
#'@export
vc_test_perm <- function(y, x, indiv=rep(1,nrow(x)), phi, w, Sigma_xi = diag(ncol(phi)),
                         n_perm=1000){

  score_obs <- vc_score(y = y, x = x, indiv = indiv, phi = phi, w = w, Sigma_xi = Sigma_xi)$score
  n_samples <- ncol(y)

  scores_perm <- numeric(n_perm)
  indiv_fact <- factor(indiv)

  if(is.null(colnames(y))){
    colnames(y) <- 1:ncol(y)
  }

  strat_sampling <- function(fact){
    res <- numeric(length(fact))
    for(l in levels(fact)){
      original_index <- which(fact==l)
      res[original_index] <- as.numeric(sample(as.character(original_index)))
    }
    return(res)
  }

  for(b in 1:n_perm){
    ## permute samples within indiv
    perm_index <- strat_sampling(indiv_fact)
    scores_perm[b] <- vc_score(y[, perm_index], x[perm_index, , drop=FALSE], indiv_fact, phi, w, Sigma_xi = Sigma_xi)$score
  }

  pval <- 1-sum(scores_perm < score_obs)/n_perm

  return(list("scores_perm" = scores_perm, "score_obs" = score_obs, "pval" = pval))

}
