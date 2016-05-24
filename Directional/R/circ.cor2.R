################################
#### Cicrular correlation for cicular data II
#### Tsagris Michail 3/2014
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics
################################

circ.cor2 <- function(theta, phi, rads = FALSE) {
  ## theta and phi are angular data in degrees or radians
  ## by default they are in degrees
  n <- length(theta)  ## sample size
  ## if the data are in degrees we transform them into radians
  if (rads == FALSE) {
    theta <- theta * pi/180
    phi <- phi * pi/180
  }
  rcc <- cor( cos(theta), cos(phi) )
  rcs <- cor( cos(theta), sin(phi) )
  rss <- cor( sin(theta), sin(phi) )
  rsc <- cor( sin(theta), cos(phi) )
  r1 <- cor( cos(theta), sin(theta) )
  r2 <- cor( cos(phi), sin(phi) )
  up <- rcc^2 + rcs^2 + rsc^2 + rss^2 + 2 * (rcc * rss + rcs * rsc) * r1 * r2 -
    2 * (rcc * rcs + rsc * rss) * r2 - 2 * (rcc * rsc + rcs * rss) * r1
  down <- (1 - r1^2) * (1 - r2^2)
  rho <- up/down
  test <- n * rho^2
  pvalue <- 1 - pchisq(test, 4)
  res <- c(rho, pvalue)
  names(res) <- c("rho", "p-value")
  res
}
