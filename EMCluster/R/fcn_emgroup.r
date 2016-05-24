### This file contains all wraps to call C function in "src/R_M_emgroup.c".
### Written: Wei-Chen Chen on 2008/10/07.


# Call:
#   SEXP R_M_emgroup(SEXP x, SEXP n, SEXP p, SEXP nclass,
#                    SEXP alpha, SEXP em_iter, SEXP em_eps)
# Input:
#   x: SEXP[n, p], data matrix of n*p.
#   n: SEXP[1], number of observations.
#   p: SEXP[1], number of dimersions.
#   nclass: SEXP[1], number of classes.
#   alpha: SEXP[1], 0.99 by default.
#   em_iter: SEXP[1], number of iterations, 1000 by default.	# em.iter
#   em_eps: SEXP[1], tolerance, 1e-4 by default.		# em.eps
# Output in C:
#   ret: a list contains
#     pi: SEXP[nclass], proportions of classes.
#     Mu: SEXP[nclass, p], means of MVNs.
#     LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular Sigma matrices.
#     llhdval: SEXP[1], log likelihood value.
#     nc: SEXP[nclass], number of observations in each class.
#     class: SEXP[n], class id's for all observations
#            starting from 0 to (nclass - 1).
#     flag: SEXP[1], returned value from emgroup() in "src/emgroup.c".
# Output in R:
#     n: SEXP[1], number of observations.
#     p: SEXP[1], number of dimersions.
#     nclass: SEXP[1], number of classes.
#     method: SEXP[1], initialization method.
emgroup <- function(x, nclass = 1, EMC = .EMC){
  n <- nrow(x)
  p <- ncol(x)

  ret <- .Call("R_M_emgroup",
               as.double(t(x)),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.double(EMC$alpha),
               as.integer(EMC$em.iter),
               as.double(EMC$em.eps))

  if(ret$flag == 1) ret$llhdval <- NA
  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, nrow = nclass, byrow = TRUE)
  ret$LTSigma <- matrix(ret$LTSigma, nrow = nclass, byrow = TRUE)
  ret$class <- ret$class + 1

  ret$n <- n
  ret$p <- p 
  ret$nclass <- nclass
  ret$method <- "svd.kmeans.em"

  class(ret) <- "emret"
  ret
}

### summary() and print.summary() for class "emret" are in "R/fcn_summary.r".




### This function is only for advance users.
emgroup.wt <- function(x, nclass = 1, EMC = .EMC){
  n <- ncol(x)
  p <- nrow(x)

  ret <- .Call("R_M_emgroup",
               as.double(x),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.double(EMC$alpha),
               as.integer(EMC$em.iter),
               as.double(EMC$em.eps))

  if(ret$flag == 1) ret$llhdval <- NA
  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, ncol = nclass)
  ret$LTSigma <- matrix(ret$LTSigma, ncol = nclass)
  ret$class <- ret$class + 1

  ret$n <- n
  ret$p <- p 
  ret$nclass <- nclass
  ret$method <- "svd.kmeans.em"

  class(ret) <- "emret.wt"
  ret
}

