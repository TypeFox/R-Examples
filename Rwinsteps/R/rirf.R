rirf <- function(x, theta = seq(-4, 4, length = 100)){

  if(class(x) == "ifile")
    x <- x$measure

  ni <- length(x)
  out <- list(theta = theta)
  out$p <- sapply(1:ni, function(i) 1/(1 + exp(-theta + x[i])))

  class(out) <- "rirf"

  return(out)
}
