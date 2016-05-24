BRInitialValue2ndLevelMeanKnown <- function(given) {
  # This function makes the initial values needed to run BRIMM.
  # "Kn" means the descriptive second level mean (mean of Beta distribution) is known.

  if (all(given$prior.mean == 0)) {
    if(all(given$sample.mean == mean(given$sample.mean))) {
      r.ini <- mean(given$sample.mean) * (1 - mean(given$sample.mean)) / (var(given$sample.mean) + 0.1)
    } else {
      r.ini <- mean(given$sample.mean) * (1 - mean(given$sample.mean)) / var(given$sample.mean)
    }
  } else {
    if(all(given$sample.mean == mean(given$sample.mean))) {
      r.ini <- mean(given$prior.mean) * (1 - mean(given$prior.mean)) / (var(given$sample.mean) + 0.1)
    } else {
      r.ini <- mean(given$prior.mean) * (1 - mean(given$prior.mean)) / var(given$sample.mean)
    }
  }
  list(r.ini = r.ini, a.ini = -log(r.ini))
}

BRInitialValue2ndLevelMeanUnknown <- function(given) {
  # This function makes the initial values needed to run BRIMM.
  # "Un" means the descriptive second level mean (mean of Beta distribution) is unknown.

  y <- given$sample.mean
  z <- given$z
  n <- given$n
  x.ini <- given$x.ini

  if (identical(x.ini, NA) & !given$intercept) {
    print("In case we do not designate prior mean, at least one beta should be estimated 
           (either intercept=T or at least one covariate is needed)")
    stop()
  } else if (identical(x.ini, NA) & given$intercept) {
    x <- matrix(1, length(y), 1)
    b.ini <- as.vector(log(mean(y) / (1 - mean(y))))
  } else if (!identical(x.ini, NA) & given$intercept) {
    x.ini <- as.matrix(x.ini)
    xname <- paste("x.ini[, ", 1 : ncol(x.ini), "]", sep = "")
    formula <- as.formula(paste("cbind(z, n - z) ~ ", paste(xname, collapse = "+")))
    b.ini <- glm(formula, family = binomial)$coefficients
    x <- as.matrix(cbind(rep(1, length(y)), x.ini))
  } else if (!identical(x.ini, NA) & !given$intercept) {
    x.ini <- as.matrix(x.ini)
    xname <- paste("x.ini[, ", 1 : ncol(x.ini), "]", sep = "")
    formula <- as.formula(paste("cbind(z, n - z) ~ ", paste(xname, collapse = "+"), "- 1"))
    b.ini <- glm(formula, family = binomial)$coefficients
    x <- x.ini
  }	

  p0.ini <- mean(exp(x %*% b.ini) / (1 + exp(x %*% b.ini)))

  if(all(y == mean(y))) {
    r.ini <- p0.ini * (1 - p0.ini) / (var(y) + 0.1)
  } else {
    r.ini <- p0.ini * (1 - p0.ini) / var(y)
  }
  list(x = x, b.ini = b.ini, a.ini = -log(r.ini))
}

BRAlphaEst2ndLevelMeanKnown <- function(given, ini) {
  # Alpha estimation of BRIMM when the second level mean is known

  z <- given$z
  n <- given$n
  k <- length(n)
  a.ini <- ini$a.ini
  p0 <- given$prior.mean
  q0 <- 1 - p0 

  BRDerivAlpha <- function(a) {
    # The first and second order derivatives of log likelihood with respect to alpha
    zap <- z + exp(-a) * p0 
    nzaq <- n - z + exp(-a) * q0
    ap <- exp(-a) * p0 
    aq <- exp(-a) * q0
    al <- exp(-a)

    const1 <- ((digamma(zap) - digamma(ap)) * p0 + (digamma(nzaq) - digamma(aq)) * q0
               + digamma(al) - digamma(n + al))
    const3 <- ((trigamma(zap) - trigamma(ap)) * p0^2 + (trigamma(nzaq) - trigamma(aq)) * q0^2
               + trigamma(al) - trigamma(n + al))

    out <- c(1 - exp(-a) * sum(const1), exp(-a) * (sum(const1) + exp(-a) * sum(const3)))
    out
  }

  dif <- 1
  eps <- 0.0001
  while (abs(dif) > eps) { 
    out1 <- BRDerivAlpha(a.ini)
    score <- out1[1]
    hessian <- out1[2]
    dif <- score / hessian
    a.ini <- a.ini - dif
  }

  list(a.new = a.ini, a.var = - 1 / hessian)
}


BRAlphaBetaEst2ndLevelMeanUnknown <- function(given, ini) {
  # Alpha and Beta estimation of BRIMM when the second level mean is unknown

  z <- given$z
  n <- given$n
  x <- ini$x
  k <- length(n)
  a.ini <- ini$a.ini
  b.ini <- ini$b.ini
  m <- ncol(x)

  BRLogLikUn <- function(a, b) {
    # Log likelihood function of alpha and beta (regression coefficients) for BRIMM 
    # when the descriptive second level mean is unknown.

    p0.hat <- exp(x %*% as.matrix(b)) / (1 + exp(x %*% as.matrix(b)))
    t0 <- p0.hat * exp(-a)
    t1 <- (1 - p0.hat) * exp(-a)
    t2 <- exp(-a)
    if (any(c(t0, t1, t2, n - z + t1, t0 + z) <= 0)) {
      print("The components of lgamma should be positive")
      stop()
    } else {
      sum(lgamma(t0 + z) - lgamma(t0) + lgamma(n - z + t1) - lgamma(t1) + lgamma(t2) - lgamma(n + t2))
    }
  }

  BRDerivBeta <- function(a, b) {
    # The first and second order derivatives of log likelihood with respect to beta
    p <- exp(x %*% as.matrix(b)) / (1 + exp(x %*% as.matrix(b)))
    q <- 1 - p
    zap <- z + exp(-a) * p 
    nzaq <- n - z + exp(-a) * q
    ap <- exp(-a) * p 
    aq <- exp(-a) * q
    vec <- rep(NA, k)
    diagm <- rep(NA, k)
    if (any(c(ap, aq, zap, nzaq) == 0)) {
      vec[(ap == 0 | aq == 0 | zap == 0 | nzaq == 0)] <- 0
      tmp <- !(ap == 0 | aq == 0 | zap == 0 | nzaq == 0)
      vec[tmp] <- (digamma(zap[tmp]) - digamma(ap[tmp]) - digamma(nzaq[tmp]) 
                   + digamma(aq[tmp])) * exp(-a) * p[tmp] * q[tmp]
    } else {
      vec <- (digamma(zap) - digamma(ap) - digamma(nzaq) + digamma(aq)) * exp(-a) * p * q
    }
    if (any(c(ap, aq, zap, nzaq) == 0)) {
      diagm[(ap == 0 | aq == 0 | zap == 0 | nzaq == 0)] <- 0
      tmp <- !(ap == 0 | aq == 0 | zap == 0 | nzaq == 0)
      diagm[tmp] <- ((trigamma(zap[tmp]) - trigamma(ap[tmp]) + trigamma(nzaq[tmp]) - trigamma(aq[tmp])) 
                    * exp(-a) * p[tmp] * q[tmp] 
                    + (digamma(zap[tmp]) - digamma(ap[tmp]) - digamma(nzaq[tmp]) + digamma(aq[tmp])) 
                    * (q[tmp] - p[tmp])) * exp(-a) * p[tmp] * q[tmp]
    } else {
      diagm <- ((trigamma(z + exp(-a) * p) - trigamma(exp(-a) * p) 
                + trigamma(n - z + exp(-a) * q) - trigamma(exp(-a) * q)) * exp(-a) * p * q +
               (digamma(z + exp(-a) * p) - digamma(exp(-a) *p)
                - digamma(n - z + exp(-a) * q) + digamma(exp(-a) * q)) * (q - p)) * exp(-a) * p * q
    }
    out <- cbind(t(x) %*% as.vector(vec), t(x) %*% (x * as.numeric(diagm)))
    out
  } 

  BRDerivBeta2order <- function(a, b) {
    # The second order derivatives of log likelihood with respect to beta
    p <- exp(x %*% as.matrix(b)) / (1 + exp(x %*% as.matrix(b)))
    q <- 1 - p
    zap <- z + exp(-a) * p 
    nzaq <- n - z + exp(-a) * q
    ap <- exp(-a) * p 
    aq <- exp(-a) * q
    diagm <- rep(NA, k)
    if (any(c(ap, aq, zap, nzaq) == 0)) {
      diagm[(ap == 0 | aq == 0 | zap == 0 | nzaq == 0)] <- 0
      tmp <- !(ap == 0 | aq == 0 | zap == 0 | nzaq == 0)
      diagm[tmp] <- ((trigamma(zap[tmp]) - trigamma(ap[tmp]) + trigamma(nzaq[tmp]) - trigamma(aq[tmp])) 
                    * exp(-a) * p[tmp] * q[tmp] 
                    + (digamma(zap[tmp]) - digamma(ap[tmp]) - digamma(nzaq[tmp]) + digamma(aq[tmp])) 
                    * (q[tmp] - p[tmp])) * exp(-a) * p[tmp] * q[tmp]
    } else {
      diagm <- ((trigamma(z + exp(-a) * p) - trigamma(exp(-a) * p) 
                + trigamma(n - z + exp(-a) * q) - trigamma(exp(-a) * q)) * exp(-a) * p * q +
               (digamma(z + exp(-a) * p) - digamma(exp(-a) *p)
                - digamma(n - z + exp(-a) * q) + digamma(exp(-a) * q)) * (q - p)) * exp(-a) * p * q
    }
    out <- t(x) %*% (x * as.numeric(diagm))
    out
  } 

  BetaHatSubAlpha <- function(a) {
    dif <- 1
    eps <- 0.0001
    n.iter <- 1
    while (any(abs(dif) > eps)) { 
      out <- BRDerivBeta(a, b.ini)
      score <- out[, 1]
      hessian <- out[, 2 : (m + 1)]
      dif <- -chol2inv(chol(-hessian)) %*% score
      b.ini <- b.ini - dif
      n.iter <- n.iter + 1
      if (n.iter > 50) {
        stop()  
        print("Alpha estimate does not converge in Newton-Raphson")
      }
    }
    list(beta.new = b.ini, beta.hessian = BRDerivBeta2order(a, b.ini))
  }

  MarginalPostAlpha <- function(a) {
    b.sub.a <-  BetaHatSubAlpha(a)$beta.new
    a + BRLogLikUn(a, b.sub.a) - 0.5 * log(det(-BRDerivBeta2order(a, b.sub.a)))
  }
  a.temp <- optim(a.ini, MarginalPostAlpha, control = list(fnscale = -1), method= "L-BFGS-B",
                  hessian = TRUE,  lower = -Inf, upper = Inf)
  a.new <- a.temp$par
  a.hess <- a.temp$hessian
  b.temp.result <- BetaHatSubAlpha(a.new)
  b.new <- b.temp.result$beta.new
  b.hessian <- b.temp.result$beta.hessian

  list(a.new = a.new, beta.new = b.new, 
       a.var = -1 / a.hess, beta.var = chol2inv(chol(-b.hessian)))
}


BRShrinkageEst <- function(a.res, given) {	
  # This function calculates the shrinkage-related estimates

  a.new <- a.res$a.new
  a.var <- a.res$a.var

  B.hat <- exp(-a.new) / (given$n + exp(-a.new))  # shriankge = B
  inv.info <- 1 / a.var
  var.B.hat <- (B.hat * (1 - B.hat))^2 / ((B.hat * (1 - B.hat)) + inv.info)
  a1.beta <- inv.info / (1 - B.hat)
  a0.beta <- inv.info / B.hat
  central3.B <- 2 * (1 - 2 * B.hat) * B.hat * (1 - B.hat) / 
                (a1.beta + a0.beta + 1) / (a1.beta + a0.beta + 2)
  B.hat.low <- qbeta((1 - given$confidence.lvl) / 2, a1.beta, a0.beta)
  B.hat.upp <- qbeta((1 + given$confidence.lvl) / 2, a1.beta, a0.beta)

  list(B.hat = B.hat, inv.info = inv.info, var.B.hat = var.B.hat, central3.B = central3.B)
}

BRPosteriorEst2ndLevelMeanKnown <- function(B.res, given) {
  # This function calculates posterior-related estimates
  # when the descriptive second level mean is known.

  B.hat <- B.res$B.hat
  var.B.hat <- B.res$var.B.hat
  p0 <- given$prior.mean
  central3.B <- B.res$central3.B
  y <- given$sample.mean
  n <- given$n

  p.hat <- y - B.hat * (y - p0)
  c1 <- y * (1 - y)
  c2 <- -3 * y^2 + 2 * (1 + p0) * y - p0
  c3 <- (y - p0) * (3 * y - 1 - p0)
  c4 <- (y - p0)^2
  d1 <- c1 + c3 * B.hat^2
  d2 <- c2 + 2 * c3 * B.hat - 3 * c4 * B.hat^2
  d3 <- c3 - 3 * c4 * B.hat
  d4 <- c4
  var.p.hat <- (d1 - d2 * B.hat - d3 * var.B.hat + d4 * central3.B) / n
  a1.beta.p <- (p.hat * (1 - p.hat) / var.p.hat - 1) * p.hat
  a0.beta.p <- (p.hat * (1 - p.hat) / var.p.hat - 1) * (1 - p.hat)
  p.hat.low <- qbeta((1 - given$confidence.lvl) / 2, a1.beta.p, a0.beta.p)
  p.hat.upp <- qbeta((1 + given$confidence.lvl) / 2, a1.beta.p, a0.beta.p)

  list(post.mean = p.hat, post.sd = sqrt(var.p.hat), 
       post.intv.low = p.hat.low, post.intv.upp = p.hat.upp, prior.mean = p0)
}

BRPosteriorEst2ndLevelMeanUnknown <- function(B.res, a.res, ini, given){
  # This function calculates posterior-related estimates
  # when the descriptive second level mean is unknown.

  x <- as.matrix(ini$x)
  n <- given$n
  y <- given$sample.mean
  beta.new <- a.res$beta.new
  beta.var <- a.res$beta.var
  B.hat <- B.res$B.hat

  xVx <- diag(x %*% beta.var %*% t(x))
  mu0 <- exp(x %*% beta.new + xVx / 2)
  b0 <- (1 + mu0) / (mu0 * (exp(xVx) - 1)) + 2
  b1 <- mu0 * (b0 - 1)

  k1p <- b1 / (b1 + b0)
  k2p <- as.vector(k1p * (1 - k1p) / (b1 + b0 + 1))
  k1b <- as.vector(B.hat)
  k2b <- as.vector(B.res$var.B.hat)
  k3b <- as.vector(B.res$central3.B)

  p0.hat <- k1p
  p.hat <- (1 - B.hat)*y + B.hat*p0.hat

  c1 <- as.vector(y - y^2 - (y - 3 * y^2) * B.hat^2 - 2 * y^2 * B.hat^3
                  - (B.hat^2 - 2 * B.hat^3) * k1p^2)
  c2 <- as.vector(-2 * y + 3 * y^2 + 2 * (y - 3 * y^2) * B.hat + 3 * y^2 * B.hat^2
                  - (-2 * B.hat + 3 * B.hat^2) * k1p^2)
  c3 <- as.vector(y - 3 * y^2 + 3 * y^2 * B.hat - (3 * B.hat - 1) * k1p^2)
  c4 <- as.vector(y^2 - k1p^2)
  c5 <- as.vector(4 * y * B.hat^3 - (4 * y - 1) * B.hat^2 + 2 * (B.hat^2 - 2 * B.hat^3) * k1p)
  c6 <- as.vector(2 * (4 * y - 1) * B.hat - 6 * y * B.hat^2 - 2 * y + 1 
                  + 2 * (-2 * B.hat + 3 * B.hat^2) * k1p)
  c7 <- as.vector(4 * y - 1 - 6 * y * B.hat + 2 * (3 * B.hat - 1) * k1p)
  c8 <- as.vector(-2 * y + 2 * k1p)
  c9 <- as.vector(B.hat^2 - 2 * B.hat^3)
  c10 <- as.vector(-2 * B.hat + 3 * B.hat^2)
  c11 <- as.vector(3 * B.hat - 1)

  var.p.hat <- ((c1 + c2 * k1b + c3 * k2b + c4 * k3b + c5 * k1p + c6 * k1p * k1b
                 + c7 * k1p * k2b + c8 * k1p * k3b + c9 * k2p + c10 * k1b * k2p
                 + c11 * k2b * k2p + k3b * k2p) / n 
                + k2b * k2p + k2b * k1p^2 - 2 * y * k2b * k1p + y^2 * k2b + k1b^2 * k2p)
  a1.beta.p <- (p.hat * (1 - p.hat) / var.p.hat - 1) * p.hat
  a0.beta.p <- (p.hat * (1 - p.hat) / var.p.hat - 1) * (1 - p.hat)
  p.hat.low <- qbeta((1 - given$confidence.lvl) / 2, a1.beta.p, a0.beta.p)
  p.hat.upp <- qbeta((1 + given$confidence.lvl) / 2, a1.beta.p, a0.beta.p)

  list(post.mean = p.hat, post.sd = sqrt(var.p.hat),
       post.intv.low = p.hat.low, post.intv.upp = p.hat.upp, prior.mean = p0.hat)
}


######
BRAR <- function(given, ini, n.AR = n.AR, trial.scale = trial.scale, 
                 n.AR.factor = n.AR.factor, t, u) {

  logpost <- function(para) {
    a <- para[1]
    b <- para[2 : length(para)]
    if (m == 1) {
      p0.hat <- exp(x * b) / (1 + exp(x * b)) 
    } else {
      p0.hat <- exp(x %*% as.vector(b)) / (1 + exp(x %*% as.vector(b))) 
    }
    t0 <- p0.hat * exp(-a)
    t1 <- (1 - p0.hat) * exp(-a)
    t2 <- exp(-a)
    if (any(c(t0, t1, t2, n - z + t1, t0 + z) <= 0)) {
      print("The components of lgamma should be positive")
      stop()
    } else {
       loglik <- sum(lgamma(t0 + z) - lgamma(t0) + lgamma(n - z + t1) - 
                    lgamma(t1) + lgamma(t2) - lgamma(n + t2))
      -a -(u + 1) * log(t + exp(-a)) + loglik
    }
  }

  z <- given$z
  n <- given$n
  y <- given$sample.mean
  x <- ini$x
  k <- length(n)
  a.ini <- ini$a.ini
  b.ini <- ini$b.ini
  m <- ncol(x)

  initial.tuning <- optim(as.numeric(c(a.ini, b.ini)), logpost, control = list(fnscale = -1), 
                          method = "L-BFGS-B", hessian = TRUE, lower = -Inf, upper = Inf)
  alpha.temp <- initial.tuning$par[1]
  beta.temp <- initial.tuning$par[2 : (m + 1)]
  neg.inv.hess <- chol2inv(chol(-initial.tuning$hessian))
  alpha.var <- neg.inv.hess[1, 1]
  beta.var <- neg.inv.hess[2 : (m + 1), 2 : (m + 1)]

  if (is.na(trial.scale)) {
    scale <- 1.3 * sqrt(alpha.var)
  } else {
    scale <- trial.scale * sqrt(alpha.var)
  }

  Skewedt <- function(x, a, b) {
    (1 + x / sqrt(a + b + x^2))^(a + 0.5) * (1 - x / sqrt(a + b + x^2))^(b + 0.5) 
  }


  if (k >= 10) {
    b <- 2 * log(k)
    a <- b / 2
  } else {
    b <- 2 * k
    a <- k
  }


  n.sample <- n.AR.factor * n.AR
  n.accept <- 0

  while (n.accept < n.AR / 10) { 
  # prescreening step
  # because sometimes extreme weights appear, making n.accept very small, looping forever
  # this while function prevents it

    B <- rbeta(n.sample, a, b)
    alpha.ar.temp <- sqrt(a + b) * (2 * B - 1) / 2 / sqrt(B * (1 - B))
    mode.skewt <- (a - b) * sqrt(a + b) / sqrt(2 * a + 1) / sqrt(2 * b + 1)
    mode.adj <- alpha.temp - mode.skewt
    alpha.ar <- mode.adj + scale * alpha.ar.temp 
    alpha.den <- Skewedt(alpha.ar.temp, a, b) / scale
    df <- 4
    beta.ar <- rmt(n = n.sample, mean = beta.temp, S = beta.var * (df - 2) / df, df = df)
    beta.den <- dmt(beta.ar, mean = beta.temp, S = beta.var * (df - 2) / df, df = df)

    if (m == 1) {
      ab.logpost <- sapply(1 : n.sample, function(j) { 
                       logpost(c(alpha.ar[j], beta.ar[j]))
                    })
    } else {
      ab.logpost <- sapply(1 : n.sample, function(j) { 
                      logpost(c(alpha.ar[j], beta.ar[j, ]))
                    })
    }


    weight <- exp(ab.logpost - max(ab.logpost)) / beta.den / alpha.den 
    M <- max(weight)
    U <- runif(n.sample)
    n.accept <- sum(weight / M > U)
    accept.rate <- n.accept / n.sample
    weight.index <- which(weight / M > U)
  }

  while (n.accept < n.AR) {
    n.accept.temp <- 0
    n.iter <- 1

    while (n.accept.temp < (n.AR - n.accept) / 10) { 
      # to guarantee the additional samples do not have extreme weights
      n.sample2 <- 6 * (n.AR - n.accept)
      B2 <- rbeta(n.sample2, a, b)
      alpha.ar.temp2 <- sqrt(a + b) * (2 * B2 - 1) / 2 / sqrt(B2 * (1 - B2))
      alpha.ar2 <- mode.adj + scale * alpha.ar.temp2
      alpha.den2 <- Skewedt(alpha.ar.temp2, a, b) / scale
      beta.ar2 <- rmt(n = n.sample2, mean = beta.temp, S = beta.var * (df - 2) / df, df = df)
      beta.den2 <- dmt(beta.ar2, mean = beta.temp, S = beta.var * (df - 2) / df, df = df)

      if (m == 1) {
        ab.logpost2 <- sapply(1 : n.sample2, function(j) { 
                        logpost(c(alpha.ar2[j], beta.ar2[j]))
                      })
      } else {
        ab.logpost2 <- sapply(1 : n.sample2, function(j) { 
                        logpost(c(alpha.ar2[j], beta.ar2[j, ]))
                      })
      }

      weight2 <- exp(ab.logpost2 - max(ab.logpost2)) / beta.den2 / alpha.den2 
      M <- max(weight2)
      U <- runif(n.sample2)
      n.accept.temp <- sum(weight2 / M > U)
    }

    ab.logpost <- c(ab.logpost, ab.logpost2)
    alpha.den <- c(alpha.den, alpha.den2)
    beta.den <- c(beta.den, beta.den2)
    weight <- exp(ab.logpost - max(ab.logpost)) / beta.den / alpha.den 
    M <- max(weight)
    U <- runif(n.sample + n.sample2)
    n.accept <- sum(weight / M > U)
    accept.rate <- n.accept / n.sample
    weight.index <- which(weight / M > U)
    alpha.ar <- c(alpha.ar, alpha.ar2)
    if (m == 1) {
      beta.ar <- c(beta.ar, beta.ar2)
    } else {
      beta.ar <- rbind(beta.ar, beta.ar2)
    }
    n.sample <- n.sample + n.sample2
    n.iter <- n.iter + 1
    if (n.iter > 2) {
      break()
    }
  }
 
  if (n.AR < n.accept) {
    weight.index <- weight.index[1 : n.AR]
  } else {
    n.AR <- n.accept
  }

  alpha.sample <- alpha.ar[weight.index]

  if (m == 1) {
    beta.sample <- beta.ar[weight.index]
    p0.sample <- exp(beta.sample) / (1 + exp(beta.sample))
    p.sample <- matrix(rbeta(length(z) * n.AR, z + matrix(exp(-alpha.sample) * p0.sample, 
                                                        ncol = n.AR, nrow = length(z), byrow = T), 
                          n - z + matrix(exp(-alpha.sample) * (1 - p0.sample), 
                                         ncol = n.AR, nrow = length(z), byrow = T)),
                       nrow = length(z), ncol = n.AR)
  } else {
    beta.sample <- beta.ar[weight.index, ]
    p0.sample <- exp(beta.sample %*% t(x)) / (1 + exp(beta.sample %*% t(x)))
    p.sample <- matrix(rbeta(length(z) * n.AR, z + matrix(exp(-alpha.sample) * p0.sample, 
                                                        ncol = n.AR, nrow = length(z), byrow = T), 
                          n - z + matrix(exp(-alpha.sample) * (1 - p0.sample), 
                                         ncol = n.AR, nrow = length(z), byrow = T)),
                       nrow = length(z), ncol = n.AR)
  }


  if ( all(n == n[1]) ) {
    B.matrix <- exp(-alpha.sample) / (n[1] + exp(-alpha.sample))
    post.shrinkage <- mean(B.matrix)
  } else {
    B.matrix <- sapply(1 : length(n), function(k) {
      exp(-alpha.sample) / (n[k] + exp(-alpha.sample))
    })
    post.shrinkage <- colMeans(B.matrix)
  }

  post.m <- rowMeans(p.sample)
  post.sd <- apply(p.sample, 1, sd)
  if (m == 1) {
    prior.m <- mean(p0.sample)
    beta.mean <- mean(beta.sample)
    beta.var <- var(beta.sample)
  } else {
    prior.m <- apply(p0.sample, 2, mean)
    beta.mean <- colMeans(beta.sample)
    beta.var <- apply(beta.sample, 2, var)
  }

  post.intv <- apply(p.sample, 1, quantile, probs = c((1 - given$confidence.lvl) / 2, 1 / 2 + given$confidence.lvl / 2))
  post.intv.low <- post.intv[1, ]
  post.intv.upp <- post.intv[2, ]
  alpha.mean <- mean(alpha.sample)
  alpha.var <- var(alpha.sample)

  list(weight = weight[weight.index], shrinkage = as.numeric(post.shrinkage), post.mean = as.numeric(post.m), 
       post.sd = as.numeric(post.sd), prior.mean.hat = as.numeric(prior.m),
       post.intv.low = post.intv.low, post.intv.upp = post.intv.upp, beta.new = beta.mean, beta.var = beta.var,
       a.new = alpha.mean, a.var = alpha.var, 
       alpha.sample = alpha.sample, beta.sample = beta.sample, p.sample = p.sample,
       accept.rate = accept.rate, trial.scale = scale)
}
  

BRAR2ndLevelMeanKnown <- function(given, ini, n.AR = n.AR, t, u,
                                  trial.scale = trial.scale, n.AR.factor = n.AR.factor) {

  z <- given$z
  n <- given$n
  y <- given$sample.mean
  k <- length(n)
  a.ini <- ini$a.ini
  p0 <- given$prior.mean

  logpost <- function(para) {
    a <- para
    t0 <- p0 * exp(-a)
    t1 <- (1 - p0) * exp(-a)
    t2 <- exp(-a)
    if (any(c(t0, t1, t2, n - z + t1, t0 + z) <= 0)) {
      print("The components of lgamma should be positive")
      stop()
    } else {
      loglik <- sum(lgamma(t0 + z) - lgamma(t0) + lgamma(n - z + t1) - 
                    lgamma(t1) + lgamma(t2) - lgamma(n + t2))
      -a -(u + 1) * log(t + exp(-a)) + loglik      
    }
  }

  initial.tuning <- optim(a.ini, logpost, control = list(fnscale = -1), 
                          method = "L-BFGS-B", hessian = TRUE, lower = -Inf, upper = Inf)
  alpha.temp <- initial.tuning$par
  alpha.var <- -1 / initial.tuning$hessian

  if (is.na(trial.scale)) {
    scale <- 1.3 * sqrt(alpha.var)
  } else {
    scale <- trial.scale * sqrt(alpha.var)
  }

  Skewedt <- function(x, a, b) {
    (1 + x / sqrt(a + b + x^2))^(a + 0.5) * (1 - x / sqrt(a + b + x^2))^(b + 0.5) 
  }

  if (k >= 10) {
    b <- 2 * log(k)
    a <- b / 2
  } else {
    b <- 2 * k
    a <- k
  }
  n.sample <- n.AR.factor * n.AR
  n.accept <- 0

  while (n.accept < n.AR / 10) {
  # prescreening step
  # because sometimes extreme weights appear, making n.accept very small, looping forever
  # this while function prevents it

    B <- rbeta(n.sample, a, b)
    alpha.ar.temp <- sqrt(a + b) * (2 * B - 1) / 2 / sqrt(B * (1 - B))
    mode.skewt <- (a - b) * sqrt(a + b) / sqrt(2 * a + 1) / sqrt(2 * b + 1)
    mode.adj <- alpha.temp - mode.skewt
    alpha.ar <- mode.adj + scale * alpha.ar.temp 
    alpha.den <- Skewedt(alpha.ar.temp, a, b) / scale

    ab.logpost <- sapply(1 : n.sample, function(j) { 
                         logpost(alpha.ar[j])
                         })

    weight <- exp(ab.logpost - max(ab.logpost)) / alpha.den 
    M <- max(weight)
    U <- runif(n.sample)
    n.accept <- sum(weight / M > U)
    accept.rate <- n.accept / n.sample
    weight.index <- which(weight / M > U)
  }

  while (n.accept < n.AR) {
    n.accept.temp <- 0
    n.iter <- 1

    while (n.accept.temp < (n.AR - n.accept) / 10) { 
      # to guarantee the additional samples do not have extreme weights
      n.sample2 <- 6 * (n.AR - n.accept)
      B2 <- rbeta(n.sample2, a, b)
      alpha.ar.temp2 <- sqrt(a + b) * (2 * B2 - 1) / 2 / sqrt(B2 * (1 - B2))
      alpha.ar2 <- mode.adj + scale * alpha.ar.temp2
      alpha.den2 <- Skewedt(alpha.ar.temp2, a, b) / scale

      ab.logpost2 <- sapply(1 : n.sample2, function(j) { 
                        logpost(alpha.ar2[j])
                      })

      weight2 <- exp(ab.logpost2 - max(ab.logpost2)) / alpha.den2 
      M <- max(weight2)
      U <- runif(n.sample2)
      n.accept.temp <- sum(weight2 / M > U)
    }

    ab.logpost <- c(ab.logpost, ab.logpost2)
    alpha.den <- c(alpha.den, alpha.den2)
    weight <- exp(ab.logpost - max(ab.logpost)) / alpha.den 
    M <- max(weight)
    U <- runif(n.sample + n.sample2)
    n.accept <- sum(weight / M > U)
    accept.rate <- n.accept / n.sample
    weight.index <- which(weight / M > U)
    alpha.ar <- c(alpha.ar, alpha.ar2)
    n.sample <- n.sample + n.sample2
    n.iter <- n.iter + 1
    if (n.iter > 2) {
      break()
    }
  }

  if (n.AR < n.accept) {
    weight.index <- weight.index[1 : n.AR]
  } else {
    n.AR <- n.accept
  }

  alpha.sample <- alpha.ar[weight.index]

  p.sample <- matrix(rbeta(length(z) * n.AR, z + matrix(exp(-alpha.sample) * p0, 
                                                      ncol = n.AR, nrow = length(z), byrow = T), 
                        n - z + matrix(exp(-alpha.sample) * (1 - p0), 
                                       ncol = n.AR, nrow = length(z), byrow = T)),
                  nrow = length(z), ncol = n.AR)

  if ( all(n == n[1]) ) {
    B.matrix <- exp(-alpha.sample) / (n[1] + exp(-alpha.sample))
    post.shrinkage <- mean(B.matrix)
  } else {
    B.matrix <- sapply(1 : length(n), function(k) {
      exp(-alpha.sample) / (n[k] + exp(-alpha.sample))
    })
    post.shrinkage <- colMeans(B.matrix)
  }

  post.m <- rowMeans(p.sample)
  post.sd <- apply(p.sample, 1, sd)

  post.intv <- apply(p.sample, 1, quantile, probs = c((1 - given$confidence.lvl) / 2, 1 / 2 + given$confidence.lvl / 2))
  post.intv.low <- post.intv[1, ]
  post.intv.upp <- post.intv[2, ]
  alpha.mean <- mean(alpha.sample)
  alpha.var <- var(alpha.sample)

  list(weight = weight[weight.index], shrinkage = as.numeric(post.shrinkage), post.mean = as.numeric(post.m), 
       post.sd = as.numeric(post.sd), post.intv.low = post.intv.low, post.intv.upp = post.intv.upp,
       a.new = alpha.mean, a.var = alpha.var, alpha.sample = alpha.sample, p.sample = p.sample,
       accept.rate = accept.rate, trial.scale = scale)
}
    
br <- function(z, n, X, prior.mean, intercept = TRUE, confidence.lvl = 0.95, t = 0, u = 1, 
               n.AR = 0, n.AR.factor = 4, trial.scale = NA, save.result = TRUE){

  # The main function of BRIMM

  if (missing(X)) {
    X <- NA
  }

  if (missing(prior.mean)) {
    prior.mean <- NA
  }

  given <- list(z = z, n = n, sample.mean = z / n, x.ini = X,
                prior.mean = prior.mean, intercept = intercept, confidence.lvl = confidence.lvl)

  if (any(is.na(prior.mean))) {
    ini <- BRInitialValue2ndLevelMeanUnknown(given)
  } else {
    ini <- BRInitialValue2ndLevelMeanKnown(given)
  }

  if (n.AR == 0) {

    a.res <- if (any(is.na(prior.mean))) {
               BRAlphaBetaEst2ndLevelMeanUnknown(given, ini)
             } else {
               BRAlphaEst2ndLevelMeanKnown(given, ini)
             }

    B.res <- BRShrinkageEst(a.res, given)

    if (any(is.na(prior.mean))) {
      post.res <- BRPosteriorEst2ndLevelMeanUnknown(B.res, a.res, ini, given)
    } else {
      post.res <- BRPosteriorEst2ndLevelMeanKnown(B.res,given)
    }

    output <- list(sample.mean = given$sample.mean, se = given$n, prior.mean = prior.mean,
                   shrinkage = B.res$B.hat, sd.shrinkage = sqrt(B.res$var.B.hat), 
                   post.mean = post.res$post.mean, post.sd = post.res$post.sd, 
                   prior.mean.hat = post.res$prior.mean, post.intv.low = post.res$post.intv.low, 
                   post.intv.upp = post.res$post.intv.upp, model="br", X = X, 
                   beta.new = a.res$beta.new, beta.var = a.res$beta.var,
                   intercept = intercept, a.new = a.res$a.new, a.var = a.res$a.var, confidence.lvl = confidence.lvl,
                   weight = NA, p = NA)
    output

  } else {

    if (any(is.na(prior.mean))) {
      res <- BRAR(given, ini, n.AR = n.AR, trial.scale = trial.scale, 
                  n.AR.factor = n.AR.factor, t = t, u = u)
    } else {
      res <- BRAR2ndLevelMeanKnown(given, ini, n.AR = n.AR, trial.scale = trial.scale,
                                   n.AR.factor = n.AR.factor, t = t, u = u)
    }

    if (any(is.na(prior.mean))) {
      p0.mean <- res$prior.mean.hat
    } else {
      p0.mean <- given$prior.mean
    }

    if (any(is.na(prior.mean))) {
      b.mean <- res$beta.new
      b.var <- res$beta.var
      beta <- res$beta.sample
    } else {
      b.mean <- NA
      b.var <- NA
      beta <- NA
    }

    output <- if (save.result == TRUE) {
      list(sample.mean = given$sample.mean, se = given$n, prior.mean = prior.mean,
           shrinkage = res$shrinkage, post.mean = res$post.mean, post.sd = res$post.sd, 
           prior.mean.hat = p0.mean, post.intv.low = res$post.intv.low, 
           post.intv.upp = res$post.intv.upp, model = "br", X = X, trial.scale.est = res$trial.scale,
           beta.new = b.mean, beta.var = b.var, weight = res$weight, trial.scale = trial.scale,
           intercept = intercept, a.new = res$a.new, a.var = res$a.var, confidence.lvl = confidence.lvl, p = res$p.sample,
           alpha = res$alpha.sample, beta = beta, accept.rate = res$accept.rate, n.AR.factor = n.AR.factor,
           n.AR = n.AR, t = t, u = u)
    } else {
      list(sample.mean = given$sample.mean, se = given$n, prior.mean = prior.mean,
           shrinkage = res$shrinkage, post.mean = res$post.mean, post.sd = res$post.sd, 
           prior.mean.hat = p0.mean, post.intv.low = res$post.intv.low, 
           post.intv.upp = res$post.intv.upp, model = "br", X = X, 
           beta.new = b.mean, beta.var = b.var, weight = 1, trial.scale = trial.scale,
           intercept = intercept, a.new = res$a.new, a.var = res$a.var, confidence.lvl = confidence.lvl, p = 1,
           alpha = res$alpha.sample, beta = beta, n.AR = n.AR)

    }
    output

  }
}
