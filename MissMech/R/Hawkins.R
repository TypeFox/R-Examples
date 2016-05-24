Hawkins <- function(data, spatcnt)
{
# This function performs the Hawkin's method for testing equality of
# covariances among groups. It is assumed that the data in y is ordered so
# that the cases from the same group are adjacent.
#gind is a vector indicating the end of the cases for each group. For
#example [8 23 30] means that there are three groups with cases 1-8 in one
#group, 9-23 and 24-30 in another.

# Also in the output, it gives the statistic computed for each group in a
# cell array, called A. This statistic should be tested for uniformity.
# also ni(i) gives the number of components of each [a[i]]
  y <- data
  n <- nrow(y)
  p <- ncol(y)
  g <- length(spatcnt)
  spool <- matrix(0, p, p)
  gind <- c(0, spatcnt)
  ygc <- matrix(0, n, p)
  ni <- matrix(0, g, 1)
  for(i in 1:g)
  {
      yg <- y[seq(gind[i] + 1, gind[i + 1]), ]
      ni[i] <- nrow(yg)
      spool <- spool + (ni[i] - 1) * cov(yg)
      ygmean <- apply(yg, 2, mean)
      ygc[seq(gind[i]+1, gind[i + 1]), ] <-
      yg - matrix(ygmean, ni[i], p, byrow = TRUE)
  }
  spool <- spool / (n - g)
  spool <- solve(spool)
  f <- matrix(0, n, 1)
  nu <- n - g - 1
  a <- vector("list",g)
  for(i in 1:g)
  {
    vij <- ygc [seq(gind[i] + 1, gind[i + 1]), ]
    vij <- apply(vij %*% spool * vij, 1, sum)
    vij <-  vij*ni[i]
    f[seq(gind[i] + 1, gind[i + 1])] <- ((n - g - p) * vij)/
                                          (p * ((ni[i] -1 ) * (n - g) - vij))
    a[[i]] <- 1 - pf(f[seq(gind[i] + 1, gind[i + 1])], p, (nu - p + 1))
  }
 list(fij = f, a = a, ni = ni)
} 
