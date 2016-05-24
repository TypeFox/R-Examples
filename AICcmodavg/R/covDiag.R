##covariance diagnostic for N-mixture model
##code slightly modified from Dennis et al. 2015: Biometrics 71: 237-246
##values <= 0 suggest lambda is infinite (data too sparse) and is likely
##to introduce problems during model fitting

##generic
covDiag <- function(object, ...){
  UseMethod("covDiag", object)
}

covDiag.default <- function(object, ...){
  stop("\nFunction not yet defined for this object class\n")
}



##for unmarkedFramePcount
covDiag.unmarkedFramePCount <- function(object, ...){
  yMat <- object@y
  p1 <- ct <- 0
  nvisits <- ncol(yMat)
  for(i in 1:(nvisits - 1)){
    for(j in (i+1):nvisits){
      p1 <- p1 + yMat[,i]*yMat[,j]
      ct <- ct+1
    }
  }
  cov.diag <- mean(p1)/ct-mean(yMat)^2
  out <- list("cov.diag" = cov.diag,
              "message" = ifelse(cov.diag <= 0,
                "Warning: lambda is infinite, data too sparse", NULL))
  class(out) <- "covDiag"
  return(out)
}



##pcount
covDiag.unmarkedFitPCount <- function(object, ...){
  yMat <- object@data@y
  p1 <- ct <- 0
  nvisits <- ncol(yMat)
  for(i in 1:(nvisits - 1)){
    for(j in (i+1):nvisits){
      p1 <- p1 + yMat[,i]*yMat[,j]
      ct <- ct+1
    }
  }
  cov.diag <- mean(p1)/ct-mean(yMat)^2
  out <- list("cov.diag" = cov.diag,
              "message" = ifelse(cov.diag <= 0,
                "Warning: lambda is infinite, data too sparse", NULL))
  class(out) <- "covDiag"
  return(out)
}



print.covDiag <- function(x, digits = 4, ...) {
  cat("\nCovariance diagnostic: ", round(x$cov.diag, digits), "\n")
  if(!is.null(x$message)) {
    cat(x$message, "\n")
  }
}
