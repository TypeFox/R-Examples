################################
#### ANOVA for cicular data (High concentration F test)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: S Rao Jammalamadaka and A SenGupta (2001)
#### Topics in circular statistics, pages 125-127
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 135
################################

hcf.circaov <- function(u, ina, rads = FALSE) {
  ## u contains all the circular data in radians or degrees
  ## ina is an indicator variable of each sample
  n <- length(u)  ## sample size
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  ## if the data are in degrees we transform them into radians
  if (rads == FALSE) u <- u * pi/180
  x1 <- cos(u)
  x2 <- sin(u)
  Ci <- aggregate(x1, by = list(ina), sum)
  Si <- aggregate(x2, by = list(ina), sum)
  Ri <- sqrt( Ci[, 2]^2 + Si[, 2]^2 )
  ## Ri is the resultant length of each group
  V <- sum(Ri)
  C <- sum(x1)
  S <- sum(x2)
  R <- sqrt(C^2 + S^2)  ## the resultant length based on all the data
   ## Next we stimate the common concentration parameter kappa
  kappa <- circ.summary(u, rads = TRUE, plot = FALSE)$kappa
  ## kappa is the estimated concentration parameter based on all the data
  if (kappa > 2) {
    Ft <- ( (n - g) * (V - R) )/( (g - 1) * (n - V) )
    pvalue <- 1 - pf(F, g - 1, n - g)
  } else  if (kappa < 2 & kappa > 1) {
    Ft <- (1 + 3/(8 * kappa)) * ((n - g) * (V - R))/( (g - 1) * (n - V) )
    pvalue <- 1 - pf(Ft, g - 1, n - g)
  } else {
    Ft <- NA
    pvalue <- NA
  }
  res <- c(Ft, pvalue, kappa)
  names(res) <- c('test', 'p-value', 'kappa')
  res
}
