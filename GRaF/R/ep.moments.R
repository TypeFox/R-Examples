# moments of the standard normal distribution for EP on the probit GP
ep.moments <- function (y, sigma2, mu) {
  
  # 1 + sigma^{2}
  sigma2p1 <- 1 + sigma2
  
  # get the likelihood
  z <- y * (mu / sqrt(sigma2p1))
  
  # get cdf and pdf of the likelihood
  cdf_z <- pnorm(z)
  pdf_z <- dnorm(z)
  
  # new \hat{\mu} (mean; second moment)
  muhat <- mu + (y * sigma2 * pdf_z) / cdf_z

  # new \hat{\sigma^{2}} (variance; third moment)
  sigma2hat <- sigma2 -
    (sigma2 ^ 2 * pdf_z) /
    (sigma2p1 * cdf_z) *
    (z + pdf_z / cdf_z)
  
  # log of the 0th moment (cdf)
  logM0 <- log(cdf_z)
  
  # return the three elements
  return (list(logM0 = logM0,
               muhat = muhat,
               sigma2hat = sigma2hat))
}