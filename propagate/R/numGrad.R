numGrad <- function(expr, envir = .GlobalEnv) 
{
  f0 <- eval(expr, envir)
  vars <- all.vars(expr)    
  p <- length(vars)    
  x <- sapply(vars, function(a) get(a, envir)) 
  
  eps <- 1e-04
  d <- 0.1
  r <- 4
  v <- 2
  
  zero.tol <- sqrt(.Machine$double.eps/7e-07) 
  h0 <- abs(d * x) + eps * (abs(x) < zero.tol)
  
  D <- matrix(0, length(f0), p)
  Daprox <- matrix(0, length(f0), r)
        
  for (i in 1:p) {    
    h <- h0
    for (k in 1:r) {
      x1 <- x2 <- x      
      x1 <- x1 + (i == (1:p)) * h
      f1 <- eval(expr, as.list(x1))
      x2 <- x2 - (i == (1:p)) * h
      f2 <- eval(expr, envir = as.list(x2))              
      Daprox[, k] <- (f1 - f2)/(2 * h[i])       
      h <- h/v      
    }
    for (m in 1:(r - 1)) for (k in 1:(r - m)) {
      Daprox[, k] <- (Daprox[, k + 1] * (4^m) - Daprox[, k])/(4^m - 1)     
    }
    D[, i] <- Daprox[, 1]        
  }
    
  return(D)
}
