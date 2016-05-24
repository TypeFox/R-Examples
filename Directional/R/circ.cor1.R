################################
#### Cicrular correlation for cicular data I
#### Tsagris Michail 1/2014
#### mtsagris@yahoo.gr
#### References: S Rao Jammalamadaka and A SenGupta (2001)
#### Topics in circular statistics
################################

circ.cor1 <- function(theta, phi, rads = FALSE) {
  ## theta and phi are angular data in degrees or radians
  ## by default they are in degrees
  n <- length(theta)  ## sample size
  ## if the data are in degrees we transform them into radians
  if (rads == F) {
    theta <- theta * pi/180
    phi <- phi * pi/180
  }
  ## We calculate the mean of each vector
  m1 <- circ.summary(theta, rads = TRUE, plot = FALSE)$mesos
  m2 <- circ.summary(phi, rads = TRUE, plot = FALSE)$mesos
  up <- sum(sin(theta - m1) * sin(phi - m2))
  down <- sqrt(sum(sin(theta - m1)^2) * sum(sin(phi - m2)^2))
  rho <- up/down  ## circular correlation
  lam22 <- mean((sin(theta - m1))^2 * (sin(phi - m2))^2)
  lam02 <- mean((sin(phi - m2))^2)
  lam20 <- mean((sin(theta - m1))^2)
  zrho <- sqrt(n) * sqrt(lam02 * lam20/lam22) * rho
  pvalue <- 2 * (1 - pnorm(abs(zrho)))
  res <- c(rho, pvalue)
  names(res) <- c("rho", "p-value")
  res
}
