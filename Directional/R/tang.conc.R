################################
#### ANOVA for cicular data (Tangential approach for equality of concentration parameters)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 141
################################

tang.conc <- function(u, ina, rads = F) {
  ## u contains all the circular data in radians or degrees
  ## ina is an indicator variable of each sample
  n <- length(u)  ## sample size
  ni <- as.vector(table(ina))
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  ## if the data are in degrees we transform them into radians
  if (rads == F) u <- u * pi/180 
  d <- NULL  ## will store the absolute sinus centred data here
  d2 <- dmi <- numeric(g)
  x1 <- cos(u)
  x2 <- sin(u)
  C <- aggregate(x1, by = list(ina), sum)
  S <- aggregate(x2, by = list(ina), sum)
  C <- as.vector( C[, -1] )
  S <- as.vector( S[, -1] )
  mi <- atan(S/C) + pi * as.numeric(C<0)
  for (i in 1:g) {
    b <- abs(sin(u[ina == i] - mi[i]))
    d <- c(d, b)
  }
  for (i in 1:g) {
    dmi[i] <- mean(d[ina == i])
    d2[i] <- sum((d[ina == i] - dmi[i])^2)
  }
  mdm <- mean(d)
  F <- ( (n - g) * sum(ni * (dmi - mdm)^2) ) / ( (g - 1) * sum(d2) )
  pvalue <- 1 - pf(F, g - 1, n - g)
  res <- c(F, pvalue)
  names(res) <- c('test', 'p-value')
  res
}