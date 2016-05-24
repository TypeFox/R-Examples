################################
#### Random values simulation from a normal mixture models for compositional data
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
#### References: Ryan P. Browne, Aisha ElSherbiny and
#### Paul D. McNicholas (2015)
#### R package mixture: Mixture Models for Clustering and Classification
################################


rmixcomp <- function(n, prob, mu, sigma, type = "alr") {
  ## n is the sample size
  ## p is a vector with the mixing probabilities
  ## mu is a matrix with with the mean vectors
  ## sigma is an array with the covariance matrices
  p2 <- c(0, cumsum(prob))
  p <- ncol(mu)  ## dimensionality of the data
  u <- runif(n)
  g <- nrow(mu)  ## how many clusters are there
  ina <- as.numeric(cut(u, breaks = p2))  ## the cluster of each observation
  ina <- sort(ina)
  nu <- as.vector(table(ina))  ## frequency table of each cluster
  y <- array(dim = c(n, p, g))
  for (j in 1:g) {
    y[1:nu[j], , j] <- rmvnorm( nu[j], mu[j, ], sigma[ , , j]) 
  }
  x <- y[1:nu[1], , 1]
  for (j in 2:g) {
    x <- rbind(x, y[1:nu[j], , j])
  }  
  if (type == "alr") { 
    x1 <- cbind(1, exp(x) ) 
    x <- x1 / rowSums(x1) 
  } else {
    x1 <- tcrossprod( x, helm( p + 1) )
    x2 <- exp(x1)
    x <- x2 / rowSums( x2 ) 
  }
  ## x is the simulated data
  ## data come from the first cluster, then from the second and so on
  list(id = ina, x = x)
}