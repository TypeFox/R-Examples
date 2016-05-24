################################
#### Normal mixture models for compositional data
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
#### References: Ryan P. Browne, Aisha ElSherbiny and
#### Paul D. McNicholas (2015)
#### R package mixture: Mixture Models for Clustering and Classification
################################

mix.compnorm <- function(x, g, model, type = "alr") {
  ## x is the compositional data
  ## g is the number of components to be used
  ## model is the type of model to be used
  ## type is either 'alr' or 'ilr'
  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x/rowSums(x)  ## makes sure x is compositional
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size

  if (type == "alr") {
    y <- log(x[, -p]/x[, p])
  } else {
    y0 <- log(x)
    y1 <- y0 - rowMeans( y0 )
    y <- tcrossprod( y1, helm(p) )

  }

  mod <- mixture::gpcm(y, G = g, mnames = model, start = 0, atol = 0.01)
  param <- mod$gpar
  mu <- matrix(nrow = g, ncol = length(param[[1]]$mu))
  su <- array(dim = c(length(param[[1]]$mu), length(param[[1]]$mu), g))
  for (i in 1:g) {
    mu[i, ] <- param[[i]]$mu  ## mean vector of each component
    su[, , i] <- param[[i]]$sigma  ## covariance of each component
  }
  prob <- param$pi  ## mixing probability of each component
  colnames(mu) <- colnames(su) <- colnames(y)
  t <- matrix(nrow = n, ncol = g)
  for (j in 1:g) {
    t[, j] <- -0.5 * log(det(2 * pi * su[, , j])) -
    0.5 * mahalanobis(y, mu[j, ], su[, , j])
  }
  pij <- prob * exp(t)/rowSums(prob * exp(t))
  est <- apply(pij, 1, which.max)
  list(type = type, mu = mu, su = su, prob = prob, est = est)
}
