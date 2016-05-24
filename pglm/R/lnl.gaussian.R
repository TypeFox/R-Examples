lnl.gaussian <- function(param, y, X, id, model, link, rn, other = NULL){
  if (is.null(other)) other <- "sd"
  Ti <- as.numeric(table(id)[as.character(id)])
  T <- Ti[1]
  N <- length(y)
  n <- length(unique(id))
  names.id <- as.character(unique(id))
  Xb <- apply(X, 2, tapply, id, mean)[as.character(id), ]
  yb <- tapply(y, id, mean)[as.character(id)]
  K <- ncol(X)
  beta <- param[1L:K]
  if (other == "sd"){
    sigmu <- param[K+1L]
    sigeps <- param[K+2L]
    gamma <- sigmu^2 / sigeps^2
    sig2 <- sigeps^2
  }
  else{
    gamma <- param[K+1L]
    sig2 <- param[K+2L]
  }
  eb <- as.numeric(yb) - as.numeric(crossprod(t(Xb), beta))
  e <- y - as.numeric(crossprod(t(X), beta))
  lnL <- - 1 / 2 * (log(2 * pi) + log(sig2) + 1 / Ti * log(1 + gamma * Ti) +
                    e^2 / sig2 - gamma * Ti / (1 + gamma * Ti) * eb^2 / sig2)
  
  lnL <- sum(lnL)

  gb <- e / sig2 * X - Ti * gamma / (1 + Ti * gamma) * eb * Xb /sig2
  gg <- - 1 / (2  * (1 + Ti * gamma)) + Ti * eb^2 / (2 * (1 + gamma * Ti)^2) / sig2
  gs <- - 1 / (2 * sig2) + e^2 / (2 * sig2^2) -
    gamma * Ti / (1 + gamma * Ti) * eb^2 / (2 * sig2^2)
  gradi <- cbind(gb, gg, gs)
  if (other == "sd"){
    ogradi <- gradi
    gmu <- 2 * sigmu / sigeps^2 * gg
    geps <- 2 * sigeps * gs - 2 * sigmu^2 / sigeps^3 * gg
    gradi <- cbind(gb, gmu, geps)
  }

  hbb <- -  crossprod(X) / sig2 + crossprod(sqrt(Ti * gamma / (1 + Ti * gamma)) * Xb) /sig2
  hbl <- - apply(Ti / (1 + gamma * Ti)^2 * eb * Xb / sig2, 2, sum)
  hll <- sum(Ti / (2 * (1 + gamma * Ti)^2) - Ti^2 / (1 + gamma * Ti)^3 * eb^2 / sig2)
  hbs <- apply(- e * X / sig2^2 + Ti * gamma / (1 + Ti * gamma) * eb * Xb /sig2^2, 2, sum)
  hls <- sum(- 1 / 2 * eb^2 / sig2^2 * Ti / (1 + gamma * Ti)^2)
  hss <- sum(1 / (2 * sig2^2) - e^2 / sig2^3 +
             gamma * Ti * eb^2 / ((1 + gamma * Ti) * sig2^3))
  H <- rbind( cbind(hbb, hbl, hbs), c(hbl, hll, hls), c(hbs, hls, hss))
  OH <- H
  H[1L:K, K+1L] <- H[K+1L, 1L:K] <- H[1L:K, K+1L] * 2 * sqrt(gamma/sig2)
  H[1L:K, K+2L] <- H[K+2L, 1L:K] <- -2 * gamma / sqrt(sig2) * OH[1L:K, K+1L] +
    2 * sqrt(sig2) * OH[1L:K, K+2L]
  H[K+1L, K+1L] <- 4 * gamma / sig2 * OH[K+1L, K+1L] + 2 * sum(ogradi[, K+1L])
#    H[K+1L, K+2L] <- H[K+2L, K+1L] <- 4 * sqrt(gamma) * OH[K+1L, K+2L] -
#      2 * sqrt(gamma / sig2) * sum(ogradi[, K+1L])
  
  H[K+1L, K+2L] <- H[K+2L, K+1L] <- 2 * sqrt(gamma/sig2) *
    (2 * sqrt(sig2) * OH[K+1L, K+2L] - 2 / sqrt(sig2) * sum(ogradi[, K+1L]) -
     2 * gamma / sqrt(sig2) * OH[K+1L, K+1L])
  
  H[K+2L, K+2L] <- 4 * sig2 * OH[K+2L, K+2L] + 2 * sum(ogradi[, K+2L]) -
    4 * gamma * OH[K+1L, K+2L] + 2 * gamma / sig2 * sum(ogradi[, K+1L])

  attr(lnL, "gradient") <- gradi
  attr(lnL, "hessian") <- H
  lnL
}

