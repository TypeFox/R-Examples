########################
#### Random value generation from a Dirichlet distribution
#### or an inverted Dirichlet distribution
#### Tsagris Michail 3/2015
#### mtsagris@yahoo.gr
#### References: Ng Kan W. and Tian Ge L. and Tang Mi L. (2011)
#### Dirichlet and Related Distributions: Theory, Methods and Applications, p.g. 177
################################

rdiri <- function(n, a) {
  ##  n is the sample size
  ##  a is the parameters vector
  D <- length(a)
  y <- matrix( rgamma(n * D, a, 1), ncol = D, byrow = TRUE )
  y / rowSums(y)  ## Dirichlet simulated values
}
