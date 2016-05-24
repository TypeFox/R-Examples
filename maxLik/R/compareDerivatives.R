compareDerivatives <- function(f, grad, hess=NULL, t0, eps=1e-6,
                               print=TRUE, ...) {
### t0 - initial parameter vector
##
## 1.  Initial function and grad eval
##  
  if(print)cat("-------- compare derivatives -------- \n")
  f0 <- f(t0, ...)
  attributes(f0) <- NULL
                           # keep only array data when printing
  if(is.function(grad))
      analytic <- grad(t0, ...)
  else if(is.numeric(grad))
      analytic = grad
  else
      stop("Argument 'grad' must be either gradient function or ",
           "pre-computed numeric gradient matrix")
  out <- list(t0=t0, f.t0=f0, compareGrad = list(analytic=analytic))
#  
  if(is.null(dim(analytic))) {
      if(print)cat("Note: analytic gradient is vector. ",
        "Transforming into a matrix form\n")
    if(length(f0) > 1)
        analytic <- matrix(analytic, length(analytic), 1)
# Note: we assume t0 is a simple vector -> hence gradient
# will be a column vector
    else
        analytic <- matrix(analytic, 1, length(analytic))
# f returns a scalar -> we have row vector along t0
  }
  if(print) {
     cat("Function value:\n")
     print(f0)
  }
  if(print)cat("Dim of analytic gradient:", dim(analytic), "\n")
  numeric <- numericGradient(f, t0, eps, ...)
  out$compareGrad$numeric = numeric 
  if(print)cat("       numeric          :", dim(numeric), "\n")
#  rDiff <- (analytic - numeric)/analytic
  rDiff <- ((analytic - numeric) /
            (0.5*(abs(analytic) + abs(numeric))) )
  rDiff[(analytic==0) & (numeric==0)] <- 0 
  rDiff. <- max(abs(rDiff), na.rm=TRUE)
  out$compareGrad$rel.diff <- rDiff 
  out$maxRelDiffGrad <- rDiff.
#
  if(print){
    if(ncol(analytic) < 2) {
      a <- cbind(t0, analytic, numeric, rDiff)
      dimnames(a) <- list(param=names(f0),
                          c("theta 0", "analytic", "numeric", "rel.diff"))
      print(a)
    }
    else {
      cat("t0\n")
      print(t0)
      cat("analytic gradient\n")
      print(analytic)
      cat("numeric gradient\n")
      print(numeric)      
      cat(paste("(anal-num)/(0.5*(abs(anal)+abs(num)))\n"))
      print(rDiff)
      a=list(t0=t0, analytic=analytic,
        numeric=numeric, rel.diff=rDiff) 
    }
    cat("Max relative difference:", rDiff., "\n")
  }
#  out <- list(t0=t0, f.t0=f0, compareGrad=a, maxRelDiffGrad=rDiff.) 
##
## Hessian?
##  
  if(!is.null(hess)) {
    if(print)cat("Comparing hessians: relative dfference\n")
    anHess <- hess(t0, ...)
    numHess <- numericGradient(grad, t0, eps, ...)
    rDifHess <- ((anHess-numHess) /
                 (0.5*(abs(anHess)+abs(numHess))) )
    rDifHess[(anHess==0) & (numHess==0)] <- 0 
    rDifHess. <- max(abs(rDifHess), na.rm=TRUE)       
    if(print)print(rDifHess.)
    out$compareHessian <- list(analytic = anHess,
                              numeric = numHess,
                              rel.diff = rDifHess)
    out$maxRelDiffHess = rDifHess.
  }
  if(print)cat("-------- END of compare derivatives -------- \n")
  invisible(out)
}
