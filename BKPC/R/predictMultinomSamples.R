predictMultinomSamples <-
function(x, beta, n.class = 1, n.burnin = 0){
  
  
  e1 <- new.env()
  x <- as.matrix(x)
  n.feat <- dim(x)[2]
  n <- dim(x)[1]
  beta <- as.matrix(beta)
  n.samples <- dim(beta)[1]
  
  p <- matrix(0, n * n.class, n.samples)
  v <- matrix(0, n, n.samples)
  
  
  for(i in 1 : n.class){ 
    SumBX <- x %*% t(beta[ ,(n.feat * (i - 1) + 1) : (n.feat * i),  drop = FALSE])
    v <- v + exp(SumBX)
    p[(n * (i - 1) + 1):(n * i), ] <- exp(SumBX)
  }  
  v <- 1 + v
  
  for(i in 1 : n.class)p[(n * (i - 1) + 1):(n * i),] <- p[(n * (i - 1) + 1):(n * i), ]/v

  
  if (n.burnin >= n.samples) stop("error: too many burn-in iterations specified")
  if (n.burnin < 0) n.burnin <- 0
  
  else if ((n.burnin + 1) == n.samples)super <- matrix(p, n, n.class)
  else super <- matrix(apply(p[ , (n.burnin + 1) : n.samples], 1, mean), n, n.class)
  e1$map <- cbind(1 - apply(super, 1, sum), super) 
  e1$class <- max.col(e1$map == apply(e1$map, 1, max))
  e1$p <- p[, (n.burnin + 1) : n.samples, drop = FALSE] 
  
  
  e1 <- as.list(e1)
  return(e1)
}
