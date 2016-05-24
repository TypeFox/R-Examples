### This file contains all wraps to call C function in "src/R_initials.c"
### and "src/ss_R_initials.c".
### Modified: Wei-Chen Chen on 2009/03/15.


# Call:
#   SEXP R_assign(SEXP x, SEXP n, SEXP p, SEXP nclass, SEXP p_LTSigma,
#                 SEXP pi, SEXP Mu, SEXP LTSigma)
# Input:
#   x: SEXP[n, p], data matrix of n*p.
#   n: SEXP[1], number of observations.
#   p: SEXP[1], number of dimersions.
#   nclass: SEXP[1], number of classes.
#   p_LTSigma: SEXP[1], dimersion of LTSigma, p * (p + 1) / 2.
#   pi: SEXP[nclass], proportions of classes.
#   Mu: SEXP[nclass, p], means of MVNs.
#   LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular sigma matrices.
# Output in C:
#   ret: a list contains
#     nc: SEXP[nclass], numbers of observations in each class.
#     class: SEXP[n], class id's for all observations
#            starting from 0 to (nclass - 1).
# Output in R:
#     nc: SEXP[nclass], numbers of observations in each class.
#     nclass: SEXP[1], number of classes.
#
# Call:
#   SEXP ss_R_assign(SEXP x, SEXP n, SEXP p, SEXP nclass, SEXP p_LTSigma,
#                    SEXP pi, SEXP Mu, SEXP LTSigma, SEXP lab)
# Input:
#   lab: SEXP[1], -1 for points with unknown clusters;
#                 0,..,(labK-1) for known.
assign.class <- function(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL,
    lab = NULL, return.all = TRUE){
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
    ret <- .Call("R_assign",
                 as.double(t(x)),
                 as.integer(n),
                 as.integer(p),
                 as.integer(nclass),
                 as.integer(p.LTSigma),
                 as.double(emobj$pi),
                 as.double(t(emobj$Mu)),
                 as.double(t(emobj$LTSigma)))
  } else{
    ret <- .Call("ss_R_assign",
                 as.double(t(x)),
                 as.integer(n),
                 as.integer(p),
                 as.integer(nclass),
                 as.integer(p.LTSigma),
                 as.double(emobj$pi),
                 as.double(t(emobj$Mu)),
                 as.double(t(emobj$LTSigma)),
                 as.integer(lab))
  }

  ret$class <- ret$class + 1

  if(return.all){
    emobj$nc <- ret$nc 
    emobj$class <- ret$class
    ret <- emobj
    class(ret) <- "emret"
  }

  ret
}

### summary() and print.summary() for class "emret" are in "R/fcn_summary.r".




### This function is only for advance users.
assign.class.wt <- function(x, emobj, lab = NULL, return.all = TRUE){
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
    ret <- .Call("R_assign",
                 as.double(x),
                 as.integer(n),
                 as.integer(p),
                 as.integer(nclass),
                 as.integer(p.LTSigma),
                 as.double(emobj$pi),
                 as.double(emobj$Mu),
                 as.double(emobj$LTSigma))
  } else{
    ret <- .Call("ss_R_assign",
                 as.double(x),
                 as.integer(n),
                 as.integer(p),
                 as.integer(nclass),
                 as.integer(p.LTSigma),
                 as.double(emobj$pi),
                 as.double(emobj$Mu),
                 as.double(emobj$LTSigma),
                 as.integer(lab))
  }

  ret$class <- ret$class + 1

  if(return.all){
    emobj$nc <- ret$nc 
    emobj$class <- ret$class
    ret <- emobj
    class(ret) <- "emret.wt"
  }

  ret
}




# Call:
#   SEXP R_meandispersion(SEXP x, SEXP n, SEXP p)
# Input:
#   x: SEXP[n, p], data matrix of n*p.
#   n: SEXP[1], number of observations.
#   p: SEXP[1], number of dimersions.
#   type: SEXP[1], 0 for original version, 1 for MLE, 2 for MME.
# Output in C:
#   ret: a list contains
#     mu: SEXP[p], means of MVNs.
#     ltsigma: SEXP[p * (p + 1) / 2], lower triangular sigma matrices.
meandispersion <- function(x, type = c("MLE", "MME", "rough")){
### x: data matrix with dimension n*p.
  n <- nrow(x)
  p <- ncol(x)
  type.index <- switch(type[1], MLE = 1, MME = 2, 0)

  ret <- .Call("R_meandispersion",
               as.double(t(x)),
               as.integer(n),
               as.integer(p),
               as.integer(type.index))

  ret$mu <- matrix(ret$mu, nrow = 1)
  ret$ltsigma <- matrix(ret$ltsigma, nrow = 1)
  ret$type <- type[1]
  ret
}
