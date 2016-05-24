print.rda <- function(x, ...){
  cat("Call:\n")
  dput(x$call)
  nonzero <- x$ngene; errors <- x$error
  #dimnames(nonzero) <- dimnames(errors) <-
  #  list(round(x$alpha, 3), round(x$delta, 3))
  tmp <- list(nonzero=nonzero, errors=errors)
  print(tmp)
}
  
