
clogitboost <- function(y, x, strata, iter, rho){
  x <- as.matrix(x)
  y <- as.numeric(y)
  fit <- fitclogitboost(y, x, strata, iter, rho)
  fit$call <- match.call()
  class(fit) <- "clogitboost"
  
  fit
}
