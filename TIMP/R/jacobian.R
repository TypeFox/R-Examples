jacobian <- function(func, nm, x, ...) {

  arglist <- list(...) 
  arglist[[nm]] <- x 
  f <- do.call(func, arglist) 
  n <- length(x) 
  eps <- sqrt(.Machine$double.eps)
  df <-matrix(NA, length(f), n)
  newarglist <- arglist 
  for (i in 1:n) {
      xx <- abs(x[i]) 
      if(xx == 0) 
	      delta <- eps 
      else 
	      delta <- eps * xx
      newx <- x
      newx[i] <- newx[i] + delta
      newarglist[[nm]] <- newx      
      df[,i] <- (do.call(func, newarglist) - f)/delta 
  }
  df 
}   
