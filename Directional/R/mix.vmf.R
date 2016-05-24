################################
#### Model based clustering using mixtures of von Mises-Fisher distributions
#### Tsagris Michail 4/2015
#### mtsagris@yahoo.gr
#### References: Kurt Hornik and  Bettina Grun (2014)
#### movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
#### http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
################################

mix.vmf <- function(x, g) {
  ## x contains the data
  ## g is the number of clusters
  x <- as.matrix(x)
  x <- x/sqrt(rowSums(x^2))
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size of the data
  lik <- NULL
  ##
  lika <- matrix(nrow = n, ncol = g)
  pij <- matrix(nrow = n, ncol = g)
  ka <- ka2 <- numeric(g)
   Apk <- function(p, k) {
     besselI(k, p/2, expon.scaled = TRUE) / besselI(k, p/2 - 1, expon.scaled = TRUE)
   }

  runtime <- proc.time()
  ## Step 1
  l <- 1
  mesa <- array(dim = c(g, p, 100))
  crit <- numeric(100)
  cl <- matrix(nrow = n, ncol = 100)
  for (vim in 1:100) {
    ini <- kmeans(x, g)  ## initially a k-means for starting values
    mesa[, , vim] <- ini$centers
    cl[, vim] <- ini$cluster
    crit[vim] <- ini$betweenss/ini$totss
  }
  epi <- which.max(crit)
  w <- as.vector(table(cl[, epi]))/n  #'# initial weights
  m1 <- mesa[, , epi]
  Rk <- sqrt(rowSums(m1^2))  ## mean resultant lengths of the initical clusters
  mat2 <- mat <- m1/Rk  ## initial mean directions
  for (j in 1:g) {
    R <- Rk[j]
    k <- numeric(4)
    i <- 1
    k[i] <- R * (p - R^2)/(1 - R^2)
    i <- 2
    k[i] <- k[i - 1] - (Apk(p, k[i - 1]) - R)/(1 - Apk(p, k[i - 1])^2 -
    (p - 1)/k[i - 1] * Apk(p, k[i - 1]))
    while (abs(k[i] - k[i - 1]) > 1e-07) {
      i <- i + 1
      k[i] <- k[i - 1] - (Apk(p, k[i - 1]) - R)/(1 - Apk(p, k[i - 1])^2 -
      (p - 1)/k[i - 1] * Apk(p, k[i - 1]))
    }
    ka[j] <- k[i] ## initial concentration parameters
    lika[, j] <- (p/2 - 1) * log(ka[j]) - 0.5 * p * log(2 * pi) -
    (log(besselI(ka[j], p/2 - 1, expon.scaled = TRUE)) + ka[j]) +
     ka[j] * (x %*% mat[j,])
  }
  lik[1] <- sum(log(rowSums(w * exp(lika))))  ## initial log-likelihood
  l <- 2
  ## Step 2
  pij <- w * exp(lika)/rowSums(w * exp(lika))  ## weights at step 2
  w <- rbind( w, colMeans(pij) ) ## weights for step 2
  for (j in 1:g) {
    m1 <- colSums(pij[, j] * x)
    mat2[j, ] <- m1/sqrt(sum(m1^2))  ## mean directions at step 2
    R <- sqrt(sum(m1^2))/sum(pij[, j])  ## mean resultant lengths at step 2
    k <- numeric(4)
    i <- 1
    k[i] <- R * (p - R^2)/(1 - R^2)
    i <- 2
    k[i] <- k[i - 1] - (Apk(p, k[i - 1]) - R)/(1 - Apk(p, k[i - 1])^2 -
    (p - 1)/k[i - 1] * Apk(p, k[i - 1]))
    while (abs(k[i] - k[i - 1]) > 1e-07) {
      i <- i + 1
      k[i] <- k[i - 1] - (Apk(p, k[i - 1]) - R)/(1 - Apk(p, k[i - 1])^2 -
      (p - 1)/k[i - 1] * Apk(p, k[i - 1]))
    }
    ka2[j] <- k[i]
    lika[, j] <- (p/2 - 1) * log(ka2[j]) - 0.5 * p * log(2 * pi) -
    (log(besselI(ka2[j], p/2 - 1, expon.scaled = T)) + ka2[j]) +
    ka2[j] * (x %*% mat2[j, ])
  }
   ka <- rbind(ka, ka2) ## concentration parameters at step 2
   mat <- abind( mat, mat2, along = 3 )
  lik[2] <- sum( log( rowSums(w[2, ] * exp(lika)) ) )  ## log-likelihood at step 2
  ## Step 3 and beyond
  while (lik[l] - lik[l - 1] > 1e-05) {
    l <- l + 1
    pij <- w[l - 1, ] * exp(lika)/rowSums( w[l - 1, ] * exp(lika) )  ## weights
    w <- rbind(w, colMeans(pij) )
    ka2 <- numeric(g)
    for (j in 1:g) {
      m1 <- colSums(pij[, j] * x)
      mat2[j, ] <- m1 / sqrt( sum(m1^2) )  ## mean directions at step l
      R <- sqrt( sum(m1^2) ) / sum(pij[, j])  ## mean resultant lengths at step l
      k <- numeric(4)
      i <- 1
      k[i] <- R * (p - R^2)/(1 - R^2)
      i <- 2
      k[i] <- k[i - 1] - (Apk(p, k[i - 1]) - R)/(1 - Apk(p, k[i - 1])^2 -
      (p - 1)/k[i - 1] * Apk(p, k[i - 1]))
      while (abs(k[i] - k[i - 1]) > 1e-07) {
        i <- i + 1
        k[i] <- k[i - 1] - (Apk(p, k[i - 1]) - R)/(1 - Apk(p, k[i - 1])^2 -
        (p - 1)/k[i - 1] * Apk(p, k[i - 1]))
      }
      ka2[j] <- k[i]
      lika[, j] <- (p/2 - 1) * log(ka2[j]) - 0.5 * p * log(2 * pi) -
      ( log(besselI(ka2[j], p/2 - 1, expon.scaled = TRUE) ) + ka2[j]) +
	ka2[j] * (x %*% mat2[j, ])
    }
    ka <- rbind(ka, ka2) ## concentration parameters at step l
    mat <- abind( mat, mat2, along = 3 )
    lik[l] <- sum(log(rowSums(w[l, ] * exp(lika))))
  }  ## log-likelihood at step l
  ta <- apply(pij, 1, which.max)  ## estimated cluster of each observation
  param <- cbind( mat[, , l], ka[l, ], table(ta)/n )
  runtime <- proc.time() - runtime
  colnames(param) <- c( paste("mu", 1:p, sep = ""), 'kappa', 'probs' )
  rownames(param) <- paste("Cluster", 1:g, sep = " ")
  list(param = param, loglik = lik[l], pred = ta, runtime = runtime)
}
