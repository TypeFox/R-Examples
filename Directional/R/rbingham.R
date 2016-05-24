################################
#### Simulating from a Bingham distribution 
#### Tsagris Michail 02/2014 
#### mtsagris@yahoo.gr
#### References: A new method to simulate the Bingham and related distributions 
#### in directional data analysis with applications
#### Kent J.T., Ganeiber A.M. and Mardia K.V. (2013)
#### http://arxiv.org/pdf/1310.8110v1.pdf 
#### and 
#### Exact Bayesian Inference for the Bingham Distribution
#### C.J. Fallaize and T. Kypraios (2014)
#### http://arxiv.org/pdf/1401.2894v1.pdf
################################

######### Simulation using any symmetric A matrix

rbingham <- function(n, A) {
  p <- ncol(A)  ## dimensionality of A
  lam <- eigen(A)$values  ## eigenvalues
  V <- eigen(A)$vectors  ## eigenvectors
  lam <- lam - lam[p]
  lam <- lam[-p]
  x <- f.rbing(n, lam)$X  ## Chris and Theo's code
  tcrossprod(x, V) ## simulated data 
}