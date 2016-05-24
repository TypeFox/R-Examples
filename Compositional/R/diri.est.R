################################
#### Dirichlet distribution parameters
#### Tsagris Michail 3/2012
#### mtsagris@yahoo.gr
################################

diri.est <- function(x, type = 'mle') {
  ## x is the compositional data
  ## type indicates how to estimate parameters
  ## type = 'mle' means the classical mle case
  ## type = 'prec' means to use the precision parameter phi
  ## type = 'ent' means to use the entropy for the estimation

  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x/rowSums(x)  ## makes sure x is compositional data
  n <- nrow(x)  ## sample size
  z <- log(x)

  ## loglik is for the 'mle' type
  loglik <- function(param, z) {
     param <- exp(param)
    -( n * lgamma( sum(param) ) - n * sum( lgamma(param) ) +
     sum( z %*% (param - 1) ) )
  }

  ## diri is for the 'prec' type
  diriphi <- function(param, z) {
   phi <- exp(param[1])
   b <- c(1, exp(param[-1]) )
   b <- b / sum(b)
   f <-  -( n * lgamma(phi) - n * sum( lgamma(phi * b) ) +
   sum( z %*% (phi * b - 1) ) )
   f
  }

  ## entro is for the 'ent' type
  entro <- function(param) {
   f <- numeric( length(param) )
   param <- exp(param)
   for ( i in 1:length(f) ) {
     f[i] <- ma[i] - digamma(param[i]) + digamma( sum(param) )
   }
   f
  }

  if (type == 'mle') {
    runtime <- proc.time()
    options(warn = -1)
    da <- nlm(loglik, colMeans(x) * 10, z = z, iterlim = 10000)
    da <- nlm(loglik, da$estimate, z = z, iterlim = 10000)
    da <- nlm(loglik, da$estimate, z = z, iterlim = 10000)
    da <- optim(da$estimate, loglik, z = z, control = list(maxit = 2000),
    hessian = TRUE)

    runtime <- proc.time() - runtime
    result <- list( loglik = -da$value, param = exp(da$par),
    std = sqrt( diag( solve(da$hessian) ) ), runtime = runtime  )
  }

  if (type == 'prec') {
    runtime <- proc.time()
    options(warn = -1)
    da <- nlm(diriphi, c(10, colMeans(x)[-1]), z = z, iterlim = 2000)
    da <- nlm(diriphi, da$estimate, z = z, iterlim = 2000)
    da <- nlm(diriphi, da$estimate, z = z, iterlim = 2000, hessian = TRUE)
    da <- optim(da$estimate, diriphi, z = z, control = list(maxit = 3000),
    hessian = TRUE)
    phi <- exp(da$par[1])
    a <- c( 1, exp(da$par[-1]) )
    a <- a/sum(a)

    runtime <- proc.time() - runtime
    result <- list( loglik = -da$value, phi = phi, a = a,
    b = phi * a, runtime = runtime )
  }

  if (type == 'ent') {
    runtime <- proc.time()
    ## this requires the BB package
    ma <- colMeans(z)
    da <- BB::BBsolve(colMeans(x) * 10, entro, control =
    list(maxit = 20000, tol = 1e-10))
    da <- BBsolve( da$par, entro, control = list(maxit = 20000, tol = 1e-10) )
    da <- BBsolve( da$par, entro, control = list(maxit = 20000, tol = 1e-10) )
    da <- BBsolve( da$par, entro, control = list(maxit = 20000, tol = 1e-10) )
    param <- exp(da$par)
    lik <- n * lgamma( sum(param) ) - n * sum( lgamma(param) ) +
    sum( z %*% (param - 1) )

    runtime <- proc.time() - runtime
    result <- list( loglik = lik, param = param, runtime = runtime )
  }

  result
}


