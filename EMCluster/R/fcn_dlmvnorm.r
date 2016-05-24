### This file contains three wraps to call C function in "src/R_dlmvnorm.c".
### Written: Wei-Chen Chen on 2008/10/28.


# Call:
#   SEXP R_lnlikelihood(SEXP x, SEXP n, SEXP p, SEXP nclass, SEXP p_LTSigma,
#                       SEXP pi, SEXP Mu, SEXP LTSigma)
# Input:
#   x: SEXP[n, p], data matrix of n*p.
#   n: SEXP[1], number of observations.
#   p: SEXP[1], number of dimersion.
#   nclass: SEXP[1], number of classes.
#   p_LTSigma: SEXP[1], dimersion of LTSigma, p * (p + 1) / 2.
#   pi: SEXP[nclass], proportions of classes.
#   Mu: SEXP[nclass, p], means of MVNs.
#   LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular Sigma matrices.
# Output in C:
#   llhdval: SEXP[1], log likelihood value.
logL <- function(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL){
### x: data matrix with dimension n*p.
### emobj (list[3]): output values contains
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

  ret <- .Call("R_lnlikelihood",
               as.double(t(x)),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.integer(p.LTSigma),
               as.double(emobj$pi),
               as.double(t(emobj$Mu)),
               as.double(t(emobj$LTSigma)))
  ret
}




### This function is only for advance users.
logL.wt <- function(x, emobj){
### x: data matrix with dimension p*n.
### emobj (list[3]): output values contains
###                  pi (array[nclass]),
###                  Mu (array[p, nclass]), and
###                  LTSigma (array[p * (p + 1) / 2, nclass]).
  n <- ncol(x)
  p <- nrow(x)
  nclass <- length(emobj$pi)
  p.LTSigma <- p * (p + 1) / 2

  ret <- .Call("R_lnlikelihood",
               as.double(x),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.integer(p.LTSigma),
               as.double(emobj$pi),
               as.double(emobj$Mu),
               as.double(emobj$LTSigma))
  ret
}




# Call:
# SEXP R_mixllhd(SEXP x, SEXP p, SEXP nclass, SEXP p_LTSigma,
#                SEXP pi, SEXP Mu, SEXP LTSigma)
# Input:
#   x: SEXP[p], data vector of length p.
#   p: SEXP[1], number of dimersions.
#   nclass: SEXP[1], number of classes.		# k
#   p_LTSigma: SEXP[1], dimersion of LTSigma, p * (p + 1) / 2.
#   pi: SEXP[nclass], proportions of classes.
#   Mu: SEXP[nclass, p], means of MVNs.
#   LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular sigma
#              matrices.
# Output:
#   dmixmvn: SEXP[1], density of mixed multivariate normal distribution.
dmixmvn <- function(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL,
                    log = FALSE){
### x: data matrix with length p or n*p.
### emobj (list[3]): output values contains
###                  pi (array[nclass]),
###                  Mu (array[nclass, p]), and
###                  LTSigma (array[nclass, p * (p + 1) / 2]).
  if(is.null(emobj)){
    emobj <- list(pi = pi, Mu = Mu, LTSigma = LTSigma)
  }

  if(is.matrix(x) || is.data.frame(x)){
    p <- ncol(x)
  } else if(is.vector(x)){
    p <- length(x)
  } else{
    stop("x should be a matrix or a vector.")
  }

  nclass <- length(emobj$pi)
  p.LTSigma <- p * (p + 1) / 2
  check.dim(emobj, p, nclass, p.LTSigma)

  my.dmixmvn <- function(x){
    .Call("R_mixllhd",
          as.double(x),
          as.integer(p),
          as.integer(nclass),
          as.integer(p.LTSigma),
          as.double(emobj$pi),
          as.double(t(emobj$Mu)),
          as.double(t(emobj$LTSigma)))
  }

  if(is.matrix(x) || is.data.frame(x)){
    ret <- apply(x, 1, my.dmixmvn)
  }
  if(is.vector(x)){
    ret <- my.dmixmvn(x) 
  }
  if(log){
    ret <- log(ret)
  }
  ret
}




### This function is only for advance users.
dmixmvn.wt <- function(x, emobj, log = FALSE){
### x: data matrix with length p or n*p.
### emobj (list[3]): output values contains
###                  pi (array[nclass]),
###                  Mu (array[nclass, p]), and
###                  LTSigma (array[nclass, p * (p + 1) / 2]).
  if(is.matrix(x) || is.data.frame(x)){
    p <- nrow(x)
  } else if(is.vector(x)){
    p <- length(x)
  } else{
    stop("x should be a matrix or a vector.")
  }

  nclass <- length(emobj$pi)
  p.LTSigma <- p * (p + 1) / 2

  my.dmixmvn <- function(x){
    .Call("R_mixllhd",
          as.double(x),
          as.integer(p),
          as.integer(nclass),
          as.integer(p.LTSigma),
          as.double(emobj$pi),
          as.double(emobj$Mu),
          as.double(emobj$LTSigma))
  }

  if(is.matrix(x) || is.data.frame(x)){
    ret <- apply(x, 2, my.dmixmvn)
  }
  if(is.vector(x)){
    ret <- my.dmixmvn(x) 
  }
  if(log){
    ret <- log(ret)
  }
  ret
}




# Call:
# SEXP R_dlmvnorm(SEXP x, SEXP p, SEXP p_LTsigma, SEXP mu, SEXP LTsigma)
# Input:
#   x: SEXP[p], data vector of length p.
#   p: SEXP[1], number of dimersions.
#   p_LTsigma: SEXP[1], dimersion of LTsigma, p * (p + 1) / 2.
#   mu: SEXP[p], means of MVN.
#   LTsigma: SEXP[p * (p + 1) / 2], lower triangular sigma matrix.
# Output:
#   dlmvn: SEXP[1], log density of multivariate normal distribution.
dlmvn <- function(x, mu, LTsigma, log = TRUE){
### x: data vector with length p.
### mu: array[p].
### LTsigma: array[nclass, p * (p + 1) / 2].
  p <- length(x)
  p.LTsigma <- p * (p + 1) / 2

  if(length(mu) != p || length(LTsigma) != p.LTsigma){
    stop("Dimensions of x, mu or LTsigma do not agree!")
  }

  ret <- .Call("R_dlmvnorm",
               as.double(x),
               as.integer(p),
               as.integer(p.LTsigma),
               as.double(mu),
               as.double(LTsigma))

  if(! log){
    ret <- exp(ret)
  }
  ret
}
dmvn <- function(x, mu, LTsigma, log = FALSE){
  dlmvn(x, mu, LTsigma, log = log)
}

