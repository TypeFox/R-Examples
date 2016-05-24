### This file contains a wrap to call C function in "src/R_mstep.c".
### Written: Wei-Chen Chen on 2008/10/19.


# Call:
#   SEXP R_mstep(SEXP x, SEXP n, SEXP p, SEXP nclass, SEXP Gamma)
# Input:
#   x: SEXP[n, p], data matrix of n*p.
#   n: SEXP[1], number of observations.
#   p: SEXP[1], number of dimersions.
#   nclass: SEXP[1], number of classes.
#   Gamma: SEXP[n, p], posteriors of classes that each observation belongs to.
# Output in C:
#   ret: a list contains
#     pi: SEXP[nclass], proportions of classes.
#     Mu: SEXP[nclass, p], means of MVNs.
#     LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular Sigma matrices.
m.step <- function(x, emobj = NULL, Gamma = NULL, assign.class = FALSE){
### x: a data matrix with dimension n*p.
### emobj (list[1]): initial values contains
###                  Gamma (array[n, nclass]).
  if(is.null(emobj)){
    emobj <- list(Gamma = Gamma)
  }

  n <- nrow(x)
  p <- ncol(x)
  nclass <- ncol(emobj$Gamma)

  if(nrow(emobj$Gamma) != n){ 
    stop("Dimensions of x or Gamma do not agree!")
  }

  ret <- .Call("R_mstep",
               as.double(t(x)),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.double(t(emobj$Gamma)))

  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, nrow = nclass, byrow = TRUE)
  ret$LTSigma <- matrix(ret$LTSigma, nrow = nclass, byrow = TRUE)

  ret$n <- n
  ret$p <- p
  ret$nclass <- nclass

  if(assign.class){
    ret <- assign.class(x, ret)
  }

  class(ret) <- "emret"
  ret
}




### This function is only for advance users.
m.step.wt <- function(x, emobj, assign.class = FALSE){
### x: a data matrix with dimension p*n.
### emobj (list[1]): initial values contains
###                  Gamma (array[nclass, n]).
  n <- ncol(x)
  p <- nrow(x)
  nclass <- nrow(emobj$Gamma)

  if(ncol(emobj$Gamma) != n){ 
    stop("Dimensions of x or Gamma do not agree!")
  }

  ret <- .Call("R_mstep",
               as.double(x),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.double(emobj$Gamma))

  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, ncol = nclass)
  ret$LTSigma <- matrix(ret$LTSigma, ncol = nclass)

  ret$n <- n
  ret$p <- p
  ret$nclass <- nclass

  if(assign.class){
    ret <- assign.class.wt(x, ret)
  }

  class(ret) <- "emret.wt"
  ret
}

