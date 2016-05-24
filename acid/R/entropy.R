entropy <-
function(x, alpha=1){ # using basic Theil index as standard, 
  x <- as.numeric(x)
  n <- length(x)
  fx<- 1/n
  if (is.null(alpha)) 
    alpha <- 1
  if (alpha == 0) { # yields mean logarithmic 
    entropy <- -sum(log(x/mean(x))*fx) 
  }
  else if (alpha == 1) 
    entropy <- sum(x/mean(x)*log(x/mean(x))*fx)
  else {
    entropy <- 1/(alpha * (alpha - 1))*sum(((x/mean(x))^alpha-1)*fx)
  }
  return(entropy)
}
