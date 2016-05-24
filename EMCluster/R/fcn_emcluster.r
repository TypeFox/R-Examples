### This file contains all wraps to call C function in "src/R_emcluster.c"
### and "src/ss_R_emcluster.c".
### Modified: Wei-Chen Chen on 2009/03/15.


# Call:
#   SEXP R_emcluster(SEXP x, SEXP n, SEXP p, SEXP nclass, SEXP p_LTSigma,
#                    SEXP pi, SEXP Mu, SEXP LTSigma, SEXP em_iter, SEXP em_eps)
# Input:
#   x: SEXP[n, p], data matrix of n*p.
#   n: SEXP[1], number of observations.
#   p: SEXP[1], number of dimersion.
#   nclass: SEXP[1], number of classes.
#   p_LTSigma: SEXP[1], dimersion of LTSigma, p * (p + 1) / 2.
#   pi: SEXP[nclass], proportions of classes.
#   Mu: SEXP[nclass, p], means of MVNs.
#   LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular Sigma matrices.
#   em.iter: SEXP[1], number of iterations, 1000 by default.	# em_iter
#   em.eps: SEXP[1], tolerance, 1e-4 by default.		# em_eps
# Output in C:
#   ret: a list contains
#     pi: SEXP[nclass], proportions of classes.
#     Mu: SEXP[nclass, p], means of MVNs.
#     LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular Sigma matrices.
#     llhdval: SEXP[1], log likelihood value.
# Output in R:
#     n: SEXP[1], number of observations.
#     p: SEXP[1], number of dimersions.
#     nclass: SEXP[1], number of classes.
#
# Call:
#   SEXP ss_R_emcluster(SEXP x, SEXP n, SEXP p, SEXP nclass, SEXP p_LTSigma,
#                       SEXP pi, SEXP Mu, SEXP LTSigma, SEXP em_iter,
#                       SEXP em_eps, SEXP lab, SEXP labK)
# Input:
#   lab: SEXP[1], -1 for points with unknown clusters;
#                 0,..,(labK-1) for known.
emcluster <- function(x, emobj = NULL, pi = NULL, Mu = NULL,
    LTSigma = NULL, lab = NULL, EMC = .EMC, assign.class = FALSE){
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

  ss <- FALSE
  if(! is.null(lab)){
    labK <- max(lab)
    lab <- lab - 1
    if(length(unique(lab[lab != -1])) != labK) stop("lab is not correct.")
    if(labK >= nclass) stop("lab is not correct.")
    if(any(table(lab[lab >= 0]) < (p + 1))) stop("lab is not correct.")
    ss <- TRUE
  }

  if(!ss){
    ret <- .Call("R_emcluster",
                 as.double(t(x)),
                 as.integer(n),
                 as.integer(p),
                 as.integer(nclass),
                 as.integer(p.LTSigma),
                 as.double(emobj$pi),
                 as.double(t(emobj$Mu)),
                 as.double(t(emobj$LTSigma)),
                 as.integer(EMC$em.iter),
                 as.double(EMC$em.eps))
  } else{
    ret <- .Call("ss_R_emcluster",
                 as.double(t(x)),
                 as.integer(n),
                 as.integer(p),
                 as.integer(nclass),
                 as.integer(p.LTSigma),
                 as.double(emobj$pi),
                 as.double(t(emobj$Mu)),
                 as.double(t(emobj$LTSigma)),
                 as.integer(EMC$em.iter),
                 as.double(EMC$em.eps),
                 as.integer(lab))
  }

  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, nrow = nclass, byrow = TRUE)
  ret$LTSigma <- matrix(ret$LTSigma, nrow = nclass, byrow = TRUE)

  ret$n <- n
  ret$p <- p
  ret$nclass <- nclass

  class(ret) <- "emret"
  if(assign.class){
    ret <- assign.class(x, ret)
    if(any(ret$nc < (p + 1))){
      warning("A stable solution is not avaliable.")
    }
  }

  ret
}

### summary() and print.summary() for class "emret" are in "R/fcn_summary.r".




### This function is only for advance users.
emcluster.wt <- function(x, emobj, lab = NULL,
   EMC = .EMC, assign.class = FALSE){
### x: data matrix with dimension p*n.
### emobj (list[3]): initial values contains
###                  pi (array[nclass]),
###                  Mu (array[p, nclass]), and
###                  LTSigma (array[p * (p + 1) / 2, nclass]).
  n <- ncol(x)
  p <- nrow(x)
  nclass <- length(emobj$pi)
  p.LTSigma <- p * (p + 1) / 2
  check.dim.wt(emobj, p, nclass, p.LTSigma)

  ss <- FALSE
  if(! is.null(lab)){
    labK <- max(lab)
    lab <- lab - 1
    if(length(unique(lab[lab != -1])) != labK) stop("lab is not correct.")
    if(labK >= nclass) stop("lab is not correct.")
    if(any(table(lab[lab >= 0]) < (p + 1))) stop("lab is not correct.")
    ss <- TRUE
  }

  if(!ss){
    ret <- .Call("R_emcluster",
                 as.double(x),
                 as.integer(n),
                 as.integer(p),
                 as.integer(nclass),
                 as.integer(p.LTSigma),
                 as.double(emobj$pi),
                 as.double(emobj$Mu),
                 as.double(emobj$LTSigma),
                 as.integer(EMC$em.iter),
                 as.double(EMC$em.eps))
  } else{
    ret <- .Call("ss_R_emcluster",
                 as.double(x),
                 as.integer(n),
                 as.integer(p),
                 as.integer(nclass),
                 as.integer(p.LTSigma),
                 as.double(emobj$pi),
                 as.double(emobj$Mu),
                 as.double(emobj$LTSigma),
                 as.integer(EMC$em.iter),
                 as.double(EMC$em.eps),
                 as.integer(lab),
                 as.integer(labK))
  }

  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, ncol = nclass)
  ret$LTSigma <- matrix(ret$LTSigma, ncol = nclass)

  ret$n <- n
  ret$p <- p
  ret$nclass <- nclass

  class(ret) <- "emret.wt"
  if(assign.class){
    ret <- assign.class.wt(x, ret)
    if(any(ret$nc < (p + 1))){
      warning("A stable solution is not avaliable.")
    }
  }

  ret
}

