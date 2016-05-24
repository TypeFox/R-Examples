#### wrapper for mvnorm
param2args <- function(param, dispstr, dim) {
  mu <- param[1:dim]
  sig2 <- param[dim + (1:dim)]
  rho <- param[-(1 : (2 * dim))]
  cormat <- diag(dim)
  if (dispstr == "ex") {
    cormat[lower.tri(cormat)] <- rho[1]
    cormat[upper.tri(cormat)] <- rho[1]
  }
  else if (dispstr == "ar1") {
    for (i in 1:dim)  for (j in 1:dim)  cormat[i,j] <- rho ^ abs(i - j)
  }
  else if (dispstr == "un") {
    cormat[lower.tri(cormat)] <- rho
    cormat[upper.tri(cormat)] <- t(cormat)[upper.tri(cormat)]
  }
  else if (dispstr == "toep") {
    for (i in 1:dim) for (j in 1:dim)
      if (i != j) cormat[i,j] <- rho[abs(i - j)]
  }
  sigmat <- diag(sqrt(sig2)) %*% cormat %*% diag(sqrt(sig2))
  list(mean = mu, sigma = sigmat, corr = cormat)
}

dmymvnorm <- function(x, param, dispstr, log=FALSE) {
  dim <- NCOL(x)
  myargs <- param2args(param, dispstr, dim)
  val <- try(dmvnorm(x, myargs$mean, myargs$sigma, log))
  if (!inherits(val, "try-error")) val else NaN
}

rmymvnorm <- function(n, param, dispstr, dim) {
  myargs <- param2args(param, dispstr, dim)
  rmvnorm(n, myargs$mean, myargs$sigma)
}
  
## pmymvnorm <- function(x, param, dispstr, seed = 123) {
##   ## seed affects dim >= 3
##   if (exists(".Random.seed", envir = .GlobalEnv)) {
##     curseed <- get(".Random.seed", .GlobalEnv)
##     on.exit(assign(".Random.seed", curseed, .GlobalEnv))
##   }
##   set.seed(seed)
##   dim <- NCOL(x)
##   myargs <- param2args(param, dispstr, dim)
##   apply(x, 1, function(x) pmvnorm(upper = x, mean = myargs$mean, sigma = myargs$sigma,
##                                   algorithm = TVPACK()))
## }

pmymvnorm <- function(x, param, dispstr) {
  dim <- NCOL(x)
  myargs <- param2args(param, dispstr, dim)
  algr <- if(dim == 2) GenzBretz() else TVPACK()
  apply(x, 1, function(x) pmvnorm(upper = x, mean = myargs$mean, sigma = myargs$sigma,
                                  algorithm = algr))
}

dmymvt <- function(x, param, dispstr, df, log=FALSE) {
  dim <- NCOL(x)
  myargs <- param2args(param, dispstr, dim)
  val <- try(dmvt(x, delta = myargs$mean, sigma = myargs$sigma, df = df, log = log, type = "shifted"))
  if (!inherits(val, "try-error")) val else NaN
}

rmymvt <- function(n, param, dispstr, df, dim) {
  myargs <- param2args(param, dispstr, dim)
  rmvt(n, sigma = myargs$sigma, df = df, delta = myargs$mean, type = "shifted")
}
  
## pmymvt <- function(x, param, dispstr, df, seed = 123) {
##   if (exists(".Random.seed", envir = .GlobalEnv)) {
##     curseed <- get(".Random.seed", .GlobalEnv)
##     on.exit(assign(".Random.seed", curseed, .GlobalEnv))
##   }
##   set.seed(seed)
##   dim <- NCOL(x)
##   myargs <- param2args(param, dispstr, dim)
##   apply(x, 1, function(x) pmvt(upper = x, df = df, delta = myargs$mean,
##                                sigma = myargs$sigma, type = "shifted",
##                                algorithm = TVPACK()))
## }

pmymvt <- function(x, param, dispstr, df) {
  dim <- NCOL(x)
  myargs <- param2args(param, dispstr, dim)
  algr <- if(dim == 2) GenzBretz() else TVPACK(1e-12)
  apply(x, 1, function(x) pmvt(upper = x, df = df, delta = myargs$mean,
                               sigma = myargs$sigma, type = "shifted",
                               algorithm = algr))
}

