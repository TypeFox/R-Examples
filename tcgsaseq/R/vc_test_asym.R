#'Computes variance component test statistic for longitudinal
#'
#'This function computes an approximation of the Variance Component test for a
#'mixture of \eqn{\chi^{2}}s using Davies method from \code{\link[CompQuadForm]{davies}}
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
#'@return A list with the following elements:\itemize{
#'   \item \code{lam}: TODO
#'   \item \code{q}: TODO
#'   \item \code{q_ext}: TODO
#'   \item \code{score_obs}: approximation of the observed score
#'   \item \code{pval}: associated p-value
#' }
#'
#'@seealso \code{\link[CompQuadForm]{davies}}
#'@importFrom stats cov
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
#'asymTestRes <- vc_test_asym(y, x, phi=t, w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
#'                            Sigma_xi=matrix(1), indiv=rep(1:4, each=3))
#'asymTestRes$pval
#'
#'@importFrom CompQuadForm davies
#'
#'@export
vc_test_asym <- function(y, x, indiv=rep(1,nrow(x)), phi, w, Sigma_xi = diag(ncol(phi))){

  score_list <- vc_score(y = y, x = x, indiv = factor(indiv), phi = phi, w = w,
                         Sigma_xi = Sigma_xi)

  Sig_q <- stats::cov(score_list$q_ext)

  if(nrow(score_list$q_ext)<2){
    warning("Only 1 individual: asymptotics likely not reached - Should probably run permutation test")
    Sig_q <- matrix(1, ncol(Sig_q), nrow(Sig_q))
  }

  lam <- svd(Sig_q)$d
  dv <- CompQuadForm::davies(score_list$score, lam)

  if(dv$ifault == 1){#error
    stop("fault in the computation from CompQuadForm::davies", dv$trace)
  }

  return(list("lam" = lam, 'q' = score_list$q, 'q_ext' = score_list$q_ext,
              "score_obs" = score_list$score, "pval" = dv$Qq)
  )
}

