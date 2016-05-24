lnl.tobit <- function(param, y, X, id, model, link, rn, other = NULL, start.sigma = FALSE){
  
  if (is.null(other)) other <- "sd"
  mills <- function(x) dnorm(x) / pnorm(x)
  if (is.null(other)) other <- "sd"
  Ti <- table(id)
  N <- length(y)
  n <- length(unique(id))
  names.id <- as.character(unique(id))
  K <- ncol(X)
  beta <- param[1L:K]
  sigma <- param[K+1L]
  if (other == "var") sigma <- sqrt(sigma)
  if (other == "lsd") sigma <- exp(sigma)
  Xb <- as.numeric(crossprod(t(X), beta))
  if (start.sigma){
    ez <- - Xb[y == 0] /sigma
    ep <- (y - Xb)[y > 0] / sigma
    mz <- mills(ez)
    fp <- fs <- numeric(length(y))
    fs[y == 0] <- - (ez + mz) * mz / sigma^2
    fs[y >  0] <- - 1 / sigma^2
    fp[y == 0] <- - mz / sigma
    fp[y >  0] <- ep / sigma
    fp <- tapply(fp, id, sum)
    fs <- tapply(fs, id, sum)
    return(sqrt(2) * sd(- fp / fs))
  }
  if (model == "pooling"){
    lnL <- numeric(length = N)
    ez <- - Xb[y == 0] /sigma
    ep <- (y - Xb)[y > 0] / sigma
    mz <- mills(ez)
    lnL[y == 0] <- log(pnorm(ez))
    lnL[y  > 0] <- - 0.5 * log(2 * pi) - log(sigma) - 0.5 *ep^2
    lnL <-  sum(lnL)
  }
  if (model == "random"){
    smu <- param[K + 2L]
    Pitr <- lapply(rn$nodes,
                   function(z){
                     result <- numeric(length = N)
                     ez <- - (Xb[y == 0] + sqrt(2) * smu * z) /sigma
                     ep <- (y - Xb - sqrt(2) * smu * z)[y > 0] / sigma
                     result[y == 0] <- pnorm(ez)
                     result[y  > 0] <- dnorm(ep) / sigma
                     result
                   }
                   )
    Pir <- lapply(Pitr, function(x) tapply(x, id, prod))
    Li <- Reduce("+", mapply("*", Pir, rn$weights, SIMPLIFY = FALSE)) / sqrt(pi)
    lnL <- sum(log(Li))
  }

  if (model == "pooling"){
    gradi <- matrix(0, nrow = nrow(X), ncol = ncol(X) + 1)
    gradi[y == 0, 1L:K] <- - mz * X[y == 0, , drop = FALSE] / sigma
    gradi[y == 0, K+1L] <- - ez * mz  / (2 * sigma^2)
    gradi[y  > 0, 1L:K] <- ep * X[y  > 0, , drop = FALSE] / sigma
    gradi[y  > 0, K+1L] <- - (1 - ep^2) / (2 * sigma^2)
    if (other == "sd") gradi[, K+1L] <- gradi[, K+1L] * (2 * sigma)
    if (other == "lsd") gradi[, K+1L] <- gradi[, K+1L] * (2 * sigma^2)
  }
  if (model == "random"){
    gradi <- Reduce("+",
                    mapply(
                           function(w, x, p){
                             ez <- - (Xb[y == 0] + sqrt(2) * smu * x) /sigma
                             ep <- (y - Xb - sqrt(2) * smu * x)[y > 0] / sigma
                             mz <- mills(ez)
                             gradi <- matrix(0, nrow = N, ncol = 2)
                             gradi[y == 0, 1] <- - mz / sigma
                             gradi[y == 0, 2] <- - ez * mz  / (2 * sigma^2)
                             gradi[y  > 0, 1] <- ep  / sigma
                             gradi[y  > 0, 2] <- - (1 - ep^2) / (2 * sigma^2)
                             gradi <- cbind(gradi, gradi[, 1] * sqrt(2) * x)
                             w * as.numeric(p[as.character(id)]) * gradi
                           },
                           rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE
                           )
                    )
    ogradi <- gradi
    if (other == "sd") gradi[, 2L] <- gradi[, 2L] * (2 * sigma)
    if (other == "lsd") gradi[, 2L] <- gradi[, 2L] * (2 * sigma^2)
    gradi <- cbind(gradi[, 1] * X, gradi[, 2:3]) /
      as.numeric(Li[as.character(id)])/ sqrt(pi)
    ogradi <- cbind(ogradi[, 1] * X, ogradi[, 2:3]) /
      as.numeric(Li[as.character(id)])/ sqrt(pi)
    
  }

  if (model == "pooling"){
    hbb <- hbs <- hss <- numeric(length = N)
    hbb[y == 0] <- - (ez + mz) * mz / sigma^2
    hbs[y == 0] <- mz * (1 - (ez + mz) * ez)/(2 * sigma^3)
    hss[y == 0] <- ez * mz * (3 - (ez + mz) * ez) / (4 * sigma^4)
    hbb[y  > 0] <- - 1 / sigma^2
    hbs[y  > 0] <- - ep / sigma^3
    hss[y  > 0] <- (1 - 2 * ep^2) / (2 * sigma^4)
    hbb <- crossprod(hbb * X, X)
    if (other == "sd"){
      hbs <- hbs * (2 * sigma)
      hss <- 4 * sigma^2 * hss + gradi[, K+1L] / sigma
    }
    if (other == "lsd"){
      hbs <- hbs * (2 * sigma^2)
      hss <- 2 * gradi[, K+1L] + 4 * sigma^4 * hss
    }
    hbs <- apply(hbs * X, 2, sum)
    hss <- sum(hss)
    H <-  rbind(cbind(hbb, hbs), c(hbs, hss))
  }
  if (model == "random"){
    H <- mapply(
                function(w, x, p){
                  P <- as.numeric((p/Li)[as.character(id)])
                  sp <- as.numeric(p / Li)
                  ez <- - (Xb[y == 0] + sqrt(2) * smu * x) /sigma
                  ep <- (y - Xb - sqrt(2) * smu * x)[y > 0] / sigma
                  mz <- mills(ez)
                  gradi <- matrix(0, nrow = N, ncol = 2)
                  gradi[y == 0, 1] <- - mz / sigma
                  gradi[y == 0, 2] <- - ez * mz  / (2 * sigma^2)
                  gradi[y  > 0, 1] <- ep  / sigma
                  gradi[y  > 0, 2] <- - (1 - ep^2) / (2 * sigma^2)
                  gradi <- cbind(gradi[, 1] * X, gradi[, 2], gradi[, 1] * sqrt(2) * x)
                  gradi <- apply(gradi, 2, tapply, id, sum)
                  H1 <- crossprod(sqrt(sp) * gradi)
                  hbb <- hbs <- hss <- numeric(length = N)
                  hbb[y == 0] <- - (ez + mz) * mz / sigma^2
                  hbs[y == 0] <- mz * (1 - (ez + mz) * ez)/(2 * sigma^3)
                  hss[y == 0] <- ez * mz * (3 - (ez + mz) * ez) / (4 * sigma^4)
                  hbb[y  > 0] <- - 1 / sigma^2
                  hbs[y  > 0] <- - ep / sigma^3
                  hss[y  > 0] <- (1 - 2 * ep^2) / (2 * sigma^4)
                  hbb <- crossprod(hbb * cbind(X, sqrt(2) * x) * P,
                                   cbind(X, sqrt(2) * x))
                  hbs <- apply(hbs * cbind(X, sqrt(2)* x) * P, 2, sum)
                  hss <- sum(hss * P)
                  H2 <- rbind(cbind(hbb, hbs), c(hbs, hss))
                  mX <- H2[1L:K, K+1L]
                  sX <- H2[1L:K, K+2L]
                  mm <- H2[K+1L, K+1L]
                  ss <- H2[K+2L, K+2L]
                  H2[K+1L, K+1L] <- ss
                  H2[K+2L, K+2L] <- mm
                  H2[1L:K, K+1L] <- H2[K+1L, 1L:K] <- sX
                  H2[1L:K, K+2L] <- H2[K+2L, 1L:K] <- mX
                  (H1 + H2) * w / sqrt(pi)
                },
                rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE
                )
    H <- Reduce("+", H) - crossprod(apply(ogradi, 2, tapply, id, sum))
    if (other == "sd"){
      H[K+1, c(1:K, K+2)] <- H[c(1:K, K+2), K+1] <- H[K+1, c(1:K, K+2)] * (2 * sigma)
      H[K+1, K+1] <- 4 * sigma^2 * H[K+1, K+1] + sum(gradi[, K+1L]) / sigma
    }
    if (other == "lsd"){
      H[K+1, c(1:K)] <- H[c(1:K), K+1] <- H[K+1, c(1:K)] * (2 * sigma^2)
      H[K+2, K+1] <- H[K+1, K+2] <- H[K+2, K+1] * (2 * sigma^2)
      H[K+1, K+1] <- 4 * sigma^4 * H[K+1, K+1] + 2 * sum(gradi[, K+1L])
    }
  }
  attr(lnL, "gradient") <- gradi
  attr(lnL, "hessian") <- H
  lnL
}
