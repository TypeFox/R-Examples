### This file contains all wraps to call C function in "src/R_init_svd.c".
### Written: Wei-Chen Chen on 2008/09/28.


# Call:
#   SEXP R_starts_via_svd(SEXP x, SEXP n, SEXP p, SEXP nclass,
#                         SEXP method, SEXP alpha)
# Input:
#   x: SEXP[n, m], a data matrix of n*m.
#   n: SEXP[1], a number of observations.
#   p: SEXP[1], a number of dimersions.
#   nclass: SEXP[1], a number of classes.
#   method: SEXP[1], an input method.
#   alpha: SEXP[1], 0.99 by default.
# Output in C:
#   ret: a list contains
#     pi: SEXP[nclass], proportions of classes.
#     Mu: SEXP[nclass, p], means of MVNs.
#     LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular sigma matrices.
#     nc: SEXP[nclass], numbers of observations in each class.
#     class: SEXP[n], class id's for all observations
#            starting from 0 to (nclass - 1).
#     flag: SEXP[1], a returned value from starts_via_svd() in
#           "src/init_svd.c".
# Output in R:
#     n: SEXP[1], a number of observations.
#     p: SEXP[1], a number of dimersions.
#     nclass: SEXP[1], a number of classes.
#     method: SEXP[1], a method.
starts.via.svd <- function(x, nclass = 1, method = c("em", "kmeans"),
    EMC = .EMC){
  if(method[1] == "em"){
    method.Call <- 1
  } else if(method[1] == "kmeans"){
    method.Call <- 0
  } else{
    stop("The method is not found.")
  }

  n <- nrow(x)
  p <- ncol(x)

  ret <- .Call("R_starts_via_svd",
               as.double(t(x)),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.integer(method.Call),
               as.double(EMC$alpha))

  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, nrow = nclass, byrow = TRUE)
  ret$LTSigma <- matrix(ret$LTSigma, nrow = nclass, byrow = TRUE)
  ret$class <- ret$class + 1

  ret$n <- n
  ret$p <- p
  ret$nclass <- nclass
  ret$method <- method

  if(ret$method == "kmeans"){
    ret$pi <- NULL
    ret$LTSigma <- NULL
  }
  class(ret) <- "svd"

  ret
}


summary.svd <- function(object, ...){
  ret <- object
  ret$mean <- t(ret$Mu)

  if(!is.null(ret$LTSigma)){
    ret$variance <- LTSigma2variance(ret$LTSigma)
  }

  ret$tp <- length(ret$pi) - 1 + length(ret$Mu) + length(ret$LTSigma)

  class(ret) <- "summary.svd"
  ret
}

print.summary.svd <- function(x,
    digits = max(4, getOption("digits") - 3), ...){
  svd <- x
  my.cat("Method: ", svd$method,
         "\n")
  my.cat(" n = ", svd$n,
         ", p = ", svd$p,
         ", nclass = ", svd$nclass,
         ", flag = ", svd$flag,
         ", total parameters = ", svd$tp,
         ".\n")
  my.cat("nc: \n")
  my.print(svd$nc, digits = digits)
  if(svd$method == "em"){
    my.cat("pi: \n")
    my.print(svd$pi, digits = digits)
  }
  my.cat("mean: \n")
  my.print(svd$mean, digits = digits)
}

print.svd <- print.summary.svd



### This function is only for advance users.
starts.via.svd.wt <- function(x, nclass = 1, method = c("em", "kmeans"),
    EMC = .EMC){
  if(method[1] == "em"){
    method.Call <- 1
  } else if(method[1] == "kmeans"){
    method.Call <- 0
  } else{
    stop("The method is not found.")
  }

  n <- ncol(x)
  p <- nrow(x)

  ret <- .Call("R_starts_via_svd",
               as.double(x),
               as.integer(n),
               as.integer(p),
               as.integer(nclass),
               as.integer(method.Call),
               as.double(EMC$alpha))

  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, ncol = nclass)
  ret$LTSigma <- matrix(ret$LTSigma, ncol = nclass)
  ret$class <- ret$class + 1

  ret$n <- n
  ret$p <- p
  ret$nclass <- nclass
  ret$method <- method

  if(ret$method == "em"){
    ret$llhdval <- logL.wt(x, ret)
    class(ret) <- "emret.wt"
  } else if(ret$method == "kmeans"){
    ret$pi <- NULL
    ret$LTSigma <- NULL
    class(ret) <- "svd.wt"
  } else{
    stop("The method is not found.")
  }
  
  ret
}

summary.svd.wt <- function(object, ...){
  ret <- wt2wot(object)
  ret$mean <- t(ret$Mu)

  if(!is.null(ret$LTSigma)){
    ret$variance <- LTSigma2variance(ret$LTSigma)
  }

  ret$tp <- length(ret$pi) - 1 + length(ret$Mu) + length(ret$LTSigma)

  class(ret) <- "summary.svd"
  ret
}
