ecmSymFit <- function(y, x, lag = 1) 
{
  if(lag < 1) stop("\nLag for ECM must be 1 or larger.\n") 
  if(!is.ts(y) || !is.ts(x) ) stop("Please provide time series data.\n")
  if (!identical(tsp(y), tsp(x))) {
    stop("Time series properties of y and x are different.\n")
  }

  name.y <- deparse(substitute(y))
  name.x <- deparse(substitute(x))
  z <- ts(data = residuals(lm(y ~ x)), start = start(y), end = end(y), 
    frequency = tsp(y)[3])  
  lz <- lag(z, k = -1) 

  dx <- diff(x); dy <- diff(y)
  xx <- bsLag(ts.union(dx, dy), lag = lag, prefix = "diff.", 
      var.name = c(name.x, name.y), suffix = ".t_", include.orig = TRUE)       
  if (tsp(xx)[1] > tsp(lz)[1]) {aa <- start(xx)} else {aa <- start(lz)} 
  if (tsp(xx)[2] < tsp(lz)[2]) {bb <- end(xx)  } else {bb <- end(lz)}   
  data <- window(cbind(xx, lz), start = aa, end = bb, frequency =tsp(y)[3])
  colnames(data) <- c(colnames(xx), "ECT.t_1")

  DepVar.x  <- data[, 1]
  DepVar.y  <- data[, lag+2]
  X.        <- data[, c(-1, -(lag+2))]
  ecm.x     <- lm(DepVar.x ~ 1 + X.)
  ecm.y     <- lm(DepVar.y ~ 1 + X.)

  result <- listn(y, x, lag, data, IndVar=X., name.y, name.x, ecm.y, ecm.x) 
  class(result) <- c("ecm", "ecmSymFit")
  return(result)
}