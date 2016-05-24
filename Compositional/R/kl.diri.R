################################
#### Kullback-Leibler divergence between two Dirichlet distributions
#### Tsagris Michail 5/2012  and
#### Bhattacharyya distance between two Dirichlet distributions
#### Tsagris Michail 7/2013  
#### mtsagris@yahoo.gr
################################

kl.diri <- function(a, b, type = "KL") {
  ## a and b are the two vectors of parameters of the two Dirichlets
  ## if type == "KL" the KL-Divergence between Dir(a) and Dir(b) is calculated
  ## if type == "bhatt" the Bhattacharyya distance between Dir(a) and 
  ## Dir(b) is calcualted
  if (type == "KL") {
    a0 <- sum(a)
    b0 <- sum(b)
    f <- sum( (a - b) * (digamma(a) - digamma(a0))) + sum(lgamma(b) - 
    lgamma(a) ) + lgamma(a0) - lgamma(b0)
  } else {
    f <- lgamma(0.5 * sum(a + b)) + 0.5 * sum(lgamma(a) + lgamma(b)) - 
	sum(lgamma(0.5 * (a + b))) - 0.5 * (lgamma(sum(a)) + lgamma(sum(b)))
  }
  f
}