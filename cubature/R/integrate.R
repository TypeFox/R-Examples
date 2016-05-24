adaptIntegrate <-
function(f, lowerLimit, upperLimit, ..., tol=1e-5,
                      fDim=1, maxEval=0, absError=0, doChecking=FALSE) {

  nL = length(lowerLimit); nU = length(upperLimit)
  if (fDim <= 0 || nL <= 0 || nU <= 0) {
      stop("Both f and x must have dimension >= 1")
  }

  if (nL != nU) {
      stop("lowerLimit and upperLimit must have same length")
  }

  if (tol <= 0) {
      stop("tol should be positive!")
  }

  f <- match.fun(f)

  if (doChecking) {
      fnF <- function(x) {
          x <- f(x, ...)
          if(!is.numeric(x) || length(x) != fDim) {
              print("adaptIntegrate: Error in evaluation function f(x) for x=")
              print(x)
              stop("adaptIntegrate: Result f(x) is not numeric or has wrong dimension")
          }
          as.double(x)
      }
  } else {
      fnF <- function(x) {
          as.double(f(x, ...))
      }
  }

  result = .Call("doCubature", as.integer(fDim), body(fnF), as.double(lowerLimit),
  as.double(upperLimit), as.integer(maxEval), as.double(absError),
  as.double(tol), environment(),
  PACKAGE="cubature")

  names(result) <- c("integral", "error", "functionEvaluations", "returnCode")
  result
}


