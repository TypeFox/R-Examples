### This file contains some tool functions.
### Written: Wei-Chen Chen on 2008/09/28.

### Functions for summary() and print.summary().
my.format <- function(x, digits = max(4, getOption("digits") - 3)){
  paste(formatC(x, format = "f", width = -1, digits = digits), collapse = " ")
}
my.cat <- function(...){
  cat(..., sep = "")
}
my.print <- function(x, digits = max(4, getOption("digits") - 3)){
  print.default(x, na.print = "", quote = FALSE, right = TRUE, digits = digits)
}


### Functions for transform objects between EMCluster and mclust.
LTsigma2var <- function(x1, p = NULL){
  if(is.null(p)){
    p <- (sqrt(1 + 8 * length(x1)) - 1) / 2
  }
  tmp <- matrix(0, nrow = p, ncol = p)
  tmp[upper.tri(tmp, diag = TRUE)] <- x1
  tmp[lower.tri(tmp)] <- t(tmp)[lower.tri(tmp)]
  tmp
}
var2LTsigma <- function(x1){
  x1[upper.tri(x1, diag = TRUE)]
}
LTSigma2variance <- function(x){
  p <- (sqrt(1 + 8 * ncol(x)) - 1) / 2
  nclass <- nrow(x)
  ret <- apply(x, 1, LTsigma2var, p)
  dim(ret) <- c(p, p, nclass)
  ret 
}
variance2LTSigma <- function(x){
  ret <- apply(x, 3, var2LTsigma)
  ret <- matrix(ret, nrow = dim(x)[3], byrow = TRUE)
  ret
}


### Function for assign Gamma by class.
class2Gamma <- function(class){
  n <- length(class)
  uniq.class <- unique(class)
  k <- length(uniq.class)

  my.z <- function(i){
    z <- vector(mode = "numeric", length = n)
    z[class == i] <- 1
    z
  }

  do.call("cbind", lapply(uniq.class, my.z))
}
class2Gamma.wt <- function(class){
  t(class2Gamma(class))
}
Gamma2class <- function(Gamma){
  unlist(apply(Gamma, 1, which.max))
}
Gamma2class.wt <- function(Gamma.wt){
  unlist(apply(Gamma.wt, 2, which.max))
}


### Function for transfer from the object class "emret.wt" to "emret".
wt2wot <- function(emobj.wt){
  ret <- emobj.wt
  if(! is.null(ret$Mu)){
    ret$Mu <- t(ret$Mu)
  }
  if(! is.null(ret$LTSigma)){
    ret$LTSigma <- t(ret$LTSigma)
  }

  if(class(emobj.wt) == "emret.wt"){
    class(ret) <- "emret"
  } else{
    class(ret) <- "emret.wt"
  }
  ret
}

check.dim <- function(emobj, p, nclass, p.LTSigma){
  if(nrow(emobj$Mu) != nclass || ncol(emobj$Mu) != p ||
     nrow(emobj$LTSigma) != nclass || ncol(emobj$LTSigma) != p.LTSigma){
    stop("Dimensions of x, pi, Mu or LTSigma do not agree!")
  }
}
check.dim.wt <- function(emobj, p, nclass, p.LTSigma){
  if(ncol(emobj$Mu) != nclass || nrow(emobj$Mu) != p ||
     ncol(emobj$LTSigma) != nclass || nrow(emobj$LTSigma) != p.LTSigma){
    stop("Dimensions of x, pi, Mu or LTSigma do not agree!")
  }
}



### Get merge results.
getlist <- function(name, alist){
  unlist(lapply(alist, function(x){ x[name] }), use.names = FALSE)
}

