fitOneDist <-
function(y, sdtype="sd.const", sd.fixed=0.05, npar)  {
  isna <- is.na(y)
  yw <- y[!isna]
  n <- length(yw)
  mu <- mean(yw)
  if (sdtype=="sd.fixed") {
    sigma <- sd.fixed
  } else {
    sigma <- sd(yw)
  }
  p <- 1.0
  psi <- list(mu=mu,sigma=sigma,p=p)
  llnw <- sum(log(sapply(yw,dnorm,psi$mu,psi$sigma)))
  NAmat <- matrix(!isna, ncol=1,  byrow=F)
  z <- matrix(nrow=length(y), ncol=1); z[NAmat] <- 1.0
  if (is.na(llnw)) llnw <- -1.0e99
  AIC <- -2*llnw + 2*npar
  BIC <- -2*llnw + log(n)*npar
  list(loglik=llnw, npar=npar, AIC=AIC, BIC=BIC, psi=psi, post=z, nobs=n, iter=1, message="")
}
