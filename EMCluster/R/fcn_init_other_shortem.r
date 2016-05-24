### This file contains a wrap to call C function in "src/R_init_other_shortem.c".
### Written: Wei-Chen Chen on 2009/01/19.


# Call:
#   SEXP R_shortemcluster(SEXP x, SEXP n, SEXP p, SEXP nclass, SEXP p_LTSigma,
#                         SEXP pi, SEXP Mu, SEXP LTSigma, SEXP maxiter,
#                         SEXP eps)
# Input:
#   x: SEXP[n, p], data matrix of n*p.
#   n: SEXP[1], number of observations.
#   p: SEXP[1], number of dimersion.
#   nclass: SEXP[1], number of classes.
#   p_LTSigma: SEXP[1], dimersion of LTSigma, p * (p + 1) / 2.
#   pi: SEXP[nclass], proportions of classes.
#   Mu: SEXP[nclass, p], means of MVNs.
#   LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular Sigma matrices.
#   maxiter: SEXP[1], number of iterations, 500 by default.	# em_iter
#   eps: SEXP[1], tolerance, 1e-2 by default.		# em_eps
# Output in C:
#   ret: a list contains
#     pi: SEXP[nclass], proportions of classes.
#     Mu: SEXP[nclass, p], means of MVNs.
#     LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular Sigma matrices.
#     llhdval: SEXP[1], log likelihood value.
#     iter: SEXP[1], iterations used in short em.
# Output in R:
#     n: SEXP[1], number of observations.
#     p: SEXP[1], number of dimersions.
#     nclass: SEXP[1], number of classes.
shortemcluster <- function(x, emobj = NULL, pi = NULL, Mu = NULL,
    LTSigma = NULL, maxiter = 100, eps = 1e-2){
### x: data matrix with dimension n*p.
### emobj (list[3]): initial values contains
###                  pi (array[nclass]),
###                  Mu (array[nclass, p]), and
###                  LTSigma (array[nclass, p * (p + 1) / 2]).
  if(is.null(emobj)){
    emobj <- list(pi = pi, Mu = Mu, LTSigma = LTSigma)
  }

  n <- nrow(x)
  p <- ncol(x)
  nclass <- length(emobj$pi)
  p.LTSigma <- p * (p + 1) / 2
  check.dim(emobj, p, nclass, p.LTSigma)
  ret <- .Call("R_shortemcluster",
               as.double(t(x)),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.integer(p.LTSigma),
               as.double(emobj$pi),
               as.double(t(emobj$Mu)),
               as.double(t(emobj$LTSigma)),
               as.integer(maxiter),
               as.double(eps))

  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, nrow = nclass, byrow = TRUE)
  ret$LTSigma <- matrix(ret$LTSigma, nrow = nclass, byrow = TRUE)

  ret$n <- n
  ret$p <- p
  ret$nclass <- nclass

  class(ret) <- "emret"
  ret
}

### summary() and print.summary() for class "emret" are in "R/fcn_summary.r".




### This function is only for advance users.
shortemcluster.wt <- function(x, emobj, maxiter = 100, eps = 1e-2){
### x: data matrix with dimension n*p.
### emobj (list[3]): initial values contains
###                  pi (array[nclass]),
###                  Mu (array[nclass, p]), and
###                  LTSigma (array[nclass, p * (p + 1) / 2]).
  n <- ncol(x)
  p <- nrow(x)
  nclass <- length(emobj$pi)
  p.LTSigma <- p * (p + 1) / 2

  ret <- .Call("R_shortemcluster",
               as.double(x),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.integer(p.LTSigma),
               as.double(emobj$pi),
               as.double(emobj$Mu),
               as.double(emobj$LTSigma),
               as.integer(maxiter),
               as.double(eps))

  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, ncol = nclass)
  ret$LTSigma <- matrix(ret$LTSigma, ncol = nclass)

  ret$n <- n
  ret$p <- p
  ret$nclass <- nclass

  class(ret) <- "emret"
  ret
}

