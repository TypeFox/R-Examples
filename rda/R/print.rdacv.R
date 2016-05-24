print.rdacv <- function(x, ...){
  cat("Call:\n")
  dput(x$call)
  nonzero <- x$ngene; errors <- x$cv.err
  #dimnames(nonzero) <- dimnames(errors) <- 
  #  list(round(x$alpha, 3), round(x$delta, 3))
  tmp <- list(nonzero=nonzero, cv.err=errors)
  print(tmp)
}
  
