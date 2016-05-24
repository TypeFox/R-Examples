################################
#### Angular central guassian
#### Tsagris Michail 2/2016
#### mtsagris@yahoo.gr
#### References: Tyler D. E. (1987). Statistical analysis for 
#### the angular central Gaussian distribution on the sphere.
#### Biometrika 74(3): 579-589.
################################

acg <- function(x) {
  p <- ncol(x)
  n <- nrow(x)
  mu <- numeric(p)
  lam <- array( dim = c(p, p, 10000) )
  lam[, , 1] <- cov(x)
  maha <- mahalanobis(x, mu, lam[, , 1])
  down <- sum( 1 / maha )
  up <- array( dim = c(p, p, n) )
  for (j in 1:n) {
    up[, , j] <- crossprod( t( x[j, ] ) ) / maha[j] 
  }
  up <- apply(up, 1:2, sum)
  lam[, , 2] <- p * up / down
  i <- 2 
  while ( sum( abs(lam[, , i] - lam[, , i - 1] ) ) > 1e-10 ) {
    i <- i + 1 
    maha <- mahalanobis(x, mu, lam[, , i - 1])
    down <- sum( 1 / maha )
    down <- sum( 1 / maha )
    up<- array(dim = c(p, p, n) )
    for (j in 1:n) {
      up[, , j] <- crossprod( t( x[j, ] ) ) / maha[j]
    }
    up <- apply(up, 1:2, sum)
    lam[, , i] <- p * up / down 
  }
  A <- lam[, , i]  
  if ( is.null( colnames(x) ) ) {
    colnames(A) <- rownames(A) <- paste("X", 1:p, sep = "")
  } else  colnames(A) <- rownames(A) <- colnames(x)
  list(iter = i, cova = A)
}

