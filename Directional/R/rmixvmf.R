
################################
#### Simulation from a model based clustering using mixtures of von Mises-Fisher distributions
#### Tsagris Michail 4/2015 
#### mtsagris@yahoo.gr
#### References: Kurt Hornik and  Bettina Grun (2014)
#### movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
#### http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
################################

rmixvmf <- function(n, prob, mu, k) {
  ## n is the sample size
  ## prob is a vector with the mixing probabilities
  ## mu is a matrix with with the mean directions
  ## k is a vector with the concentration parameters
  p2 <- c(0, cumsum(prob))
  p <- ncol(mu)  ## dimensionality of the data
  u <- runif(n)
  g <- length(k)  ## how many clusters are there
  ina <- as.numeric(cut(u, breaks = p2))  ## the cluster of each observation
  ina <- sort(ina)
  nu <- as.vector(table(ina))  ## frequency table of each cluster
  y <- array(dim = c(n, p, g))
  for (j in 1:g) {
    y[1:nu[j], , j] <- rvmf(nu[j], mu[j, ], k[j])
  }
  x <- y[1:nu[1], , 1]
  for (j in 2:g) {
    x <- rbind(x, y[1:nu[j], , j])  ## simulated data
  }  
  ## data come from the first cluster, then from the second and so on
  list(id = ina, x = x)
}