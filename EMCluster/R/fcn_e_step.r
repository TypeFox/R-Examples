### This file contains a wrap to call C function in "src/R_estep.c".
### Written: Wei-Chen Chen on 2008/10/19.


# Call:
#   SEXP R_estep(SEXP x, SEXP n, SEXP p, SEXP nclass, SEXP p_LTSigma,
#                SEXP pi, SEXP Mu, SEXP LTSigma)
# Input:
#   x: SEXP[n, p], data matrix of n*p.
#   n: SEXP[1], number of observations.
#   p: SEXP[1], number of dimersions.
#   nclass: SEXP[1], number of classes.
#   p_LTSigma: SEXP[1], dimersion of LTSigma, p * (p + 1) / 2.
#   pi: SEXP[nclass], proportions of classes.
#   Mu: SEXP[nclass, p], means of MVNs.
#   LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular Sigma matrices.
#   norm: SEXP[1], normalized.
#     TRUE returns posterior matrix, FALSE returns log(pi) + log(phi_k) matrix.
# Output in C:
#   ret: a list contains
#     Gamma: SEXP[n, p], posteriors of classes that each observation belongs to.
e.step <- function(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL,
  norm = TRUE){
### x: a data matrix with dimension n*p.
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

  ret <- .Call("R_estep",
               as.double(t(x)),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.integer(p.LTSigma),
               as.double(emobj$pi),
               as.double(t(emobj$Mu)),
               as.double(t(emobj$LTSigma)),
               as.integer(norm))

  ret$Gamma <- matrix(ret$Gamma, nrow = n, byrow = TRUE)
  ret
}




### This function is only for advance users.
e.step.wt <- function(x, emobj, norm = TRUE){
### x: a data matrix with dimension p*n.
### emobj (list[3]): initial values contains
###                  pi (array[nclass]),
###                  Mu (array[p, nclass]), and
###                  LTSigma (array[p * (p + 1) / 2, nclass]).
  n <- ncol(x)
  p <- nrow(x)
  nclass <- length(emobj$pi)
  p.LTSigma <- p * (p + 1) / 2
  check.dim.wt(emobj, p, nclass, p.LTSigma)

  ret <- .Call("R_estep",
               as.double(x),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.integer(p.LTSigma),
               as.double(emobj$pi),
               as.double(emobj$Mu),
               as.double(emobj$LTSigma),
               as.integer(norm))

  ret$Gamma <- matrix(ret$Gamma, ncol = n)
  ret
}
