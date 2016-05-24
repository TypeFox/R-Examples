### posterior samples of O-U parameters (mu, sigma, and tau)

postTheta <- function (data, X, delta, previous.theta, tau.jump, 
                        tau.thresh, tau.prior.a, tau.prior.b, sigma.prior.a, sigma.prior.b) {

  time <- data[, 1]
  leng.time <- length(time)
  time.d <- time - delta
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  time.comb <- time.temp[ord]
  time.diff <- diff(time.comb)
  leng.X <- length(X)

  mu <- previous.theta[1]
  sigma <- previous.theta[2]
  tau <- previous.theta[3]

  a.i <- exp( - time.diff / tau )   # i = 2, 3, ..., 2n
  leng.a <- length(a.i)  # 2n - 1

  # updating mu
  mu.mean <- ( X[1] + sum( (X[-1] - a.i * X[-leng.X]) / (1 + a.i) ) ) / 
                     ( 1 + sum( (1 - a.i) / (1 + a.i) ) )
  mu.sd <- sqrt( tau * sigma^2 / 2 /
                         ( 1 + sum( (1 - a.i) / (1 + a.i) ) ) )
  inv.cdf <- runif(1, min = pnorm(-30, mean = mu.mean, sd = mu.sd), 
                      max = pnorm(30, mean = mu.mean, sd = mu.sd))
  mu <- qnorm(inv.cdf, mean = mu.mean, sd = mu.sd)

  # updating sigma
  sigma <- sqrt( ( sigma.prior.b + (X[1] - mu)^2 / tau +  
                   sum( ( X[-1] - a.i * X[-leng.X] - mu * (1 - a.i) )^2 / 
                        (1 - a.i^2) ) / tau ) / 
                 rgamma(1, shape = length(time) + sigma.prior.a, scale = 1) )
 
  # updating tau
  log.post.tau <- function(t) {
    a.i.post <- exp( - time.diff / t )
    - (length(time) + 1 + tau.prior.a) * log(t) - 0.5 * sum( log(1 - a.i.post^2) ) - 
    ( tau.prior.b  + (X[1] - mu)^2 / sigma^2 +
      sum( ( X[-1] - a.i.post * X[-leng.X] - mu * (1 - a.i.post) )^2 / 
           (1 - a.i.post^2) ) / sigma^2 ) / t
  }

  tau.p <- exp(log(tau) + tau.jump)
  l.metrop <- log.post.tau(tau.p) - log.post.tau(tau)
  l.hastings <- log(tau.p) - log(tau)

  # Accept-reject
  if (l.metrop + l.hastings > tau.thresh) {
    tau <- tau.p
  }
   
  out <- c(mu, sigma, tau)  
  out
 
}



### log likelihood function of all the model parameters
logpostDelta <- function(delta, data, theta, c, log, unif, micro) {

  time <- data[, 1]
  leng.time <- length(time)

  if (delta < unif[1] | delta > unif[2]) {

    -Inf

  } else if (theta[1] < -30 | theta[1] > 30) {

    -Inf

  } else {

    lcA <- data[, 2]
    se.lcA <- data[, 3]
    lcB <- data[, 4]
    se.lcB <- data[, 5]

    if (log == TRUE) {
      # transform into magnitude scale 
      se.lcA <- se.lcA * 2.5 / lcA / log(10)
      lcA <- -2.5 * log(lcA, base = 10)
      se.lcB <- se.lcB * 2.5 / lcB / log(10)
      lcB <- -2.5 * log(lcB, base = 10)
    }

    mu <- theta[1]
    sigma <- theta[2]
    tau <- theta[3]

    # sorting time given delta
    time.d <- time - delta
    time.temp <- c(time, time.d)
    ord <- order(time.temp)
    time.comb <- time.temp[ord]
    leng.time.comb <- length(time.comb)

    # microlensing  
    if (micro == 0) {
      mat.temp <- matrix(c(rep(1, leng.time)), ncol = 1)
    } else if (micro == 1) {
      mat.temp <- matrix(c(rep(1, leng.time), time.d), ncol = 2)
    } else if (micro == 2) {
      mat.temp <- matrix(c(rep(1, leng.time), time.d, time.d^2), ncol = 3)
    } else if (micro == 3) {
      mat.temp <- matrix(c(rep(1, leng.time), time.d, time.d^2, time.d^3), ncol = 4)
    }

    c.pred <- mat.temp %*% c

    # indicator taking on 1 for X(t - delta) and 0 for X(t)
    ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]

    lc.temp <- c(lcA, lcB - c.pred)
    lc.comb <- lc.temp[ord]
    se.lc.temp <- c(se.lcA, se.lcB)
    se.lc.comb <- se.lc.temp[ord]
   
    if (all(se.lc.comb == 0)) {
      # Kelly et. al (2009)
      # x.star.i, i = 1, 2, ..., 2n
      x.star.i <- lc.comb - mu

      # omega.i, i = 1, 2, ..., 2n to be saved
      omega.i <- rep(NA, length(time.comb))

      # x.hat.i, i = 1, 2, ..., 2n to be saved
      x.hat.i <- rep(NA, length(time.comb))

      # a.i, i = 2, ..., 2n
      a.i <- exp( -diff(time.comb) / tau)

      # omega.i, i = 1, 2, ..., 2n
      omega.i[1] <- tau * sigma^2 / 2

      for (k in 2 : leng.time.comb) {
        omega.i[k] <- omega.i[1] * (1 - a.i[k - 1]^2) +
                      a.i[k - 1]^2 * omega.i[k - 1] * se.lc.comb[k - 1]^2 / 
                      (se.lc.comb[k - 1]^2 + omega.i[k - 1])
      }  

      # x.hat.i, i = 1, 2, ..., 2n
      x.hat.i[1] <- 0
      for (k in 2 : length(time.comb)) {
        x.hat.i[k] <- a.i[k - 1] * (x.hat.i[k - 1] +
                                    omega.i[k - 1] / (se.lc.comb[k - 1]^2 + omega.i[k - 1]) * 
                                    (x.star.i[k - 1] - x.hat.i[k - 1]))                  
      } 

      # log-likelihood
      sum(dnorm(x.star.i, mean = x.hat.i, sd = sqrt(omega.i + se.lc.comb^2), log = TRUE))

    } else {
      # x.star.i, i = 1, 2, ..., 2n
      x <- lc.comb - mu

      # omega.i, i = 1, 2, ..., 2n to be saved
      B <- rep(NA, leng.time.comb)

      # x.hat.i, i = 1, 2, ..., 2n to be saved
      mu.i <- rep(NA, leng.time.comb)
      mu.star.i <- rep(NA, leng.time.comb)

      # a.i, i = 2, ..., 2n
      a.i <- exp( -diff(time.comb) / tau)

      # omega.i, i = 1, 2, ..., 2n
      var0 <- tau * sigma^2 / 2
      B[1] <- se.lc.comb[1]^2 / (se.lc.comb[1]^2 + var0)

      for (k in 2 : leng.time.comb) {
        B[k] <- se.lc.comb[k]^2 / ( se.lc.comb[k]^2 + 
                                    a.i[k - 1]^2 * (1 - B[k - 1]) * se.lc.comb[k - 1]^2 +
                                    var0 * (1 - a.i[k - 1]^2) ) 
      }  
  
      # x.hat.i, i = 1, 2, ..., 2n
      mu.i[1] <- (1 - B[1]) * x[1]
      for (k in 2 : leng.time.comb) {
        mu.i[k] <- (1 - B[k]) * x[k] + B[k] * a.i[k - 1] * mu.i[k - 1]
      }
  
      mu.star.i[1] <- 0
      mu.star.i[2 : leng.time.comb] <- a.i * mu.i[-leng.time.comb]

      var.star.i <- se.lc.comb^2 / B

      # log-likelihood
      sum(dnorm(x, mean = mu.star.i, sd = sqrt(var.star.i), log = TRUE))
    }  # end of if (all(se.lc.comb == 0)) 
  }  # end of if (delta < unif[1] | delta > unif[2])
}  # end of function



### posterior samples of the latent true magnitudes
postX <- function(data, X, theta, delta, c, log, micro) {

  time <- data[, 1]
  leng.time <- length(time)

  lcA <- data[, 2]
  se.lcA <- data[, 3]
  lcB <- data[, 4]
  se.lcB <- data[, 5]

  if (log == TRUE) {
    # transform into log scale
    se.lcA <- se.lcA * 2.5 / lcA / log(10)
    lcA <- -2.5 * log(lcA, base = 10)
    se.lcB <- se.lcB * 2.5 / lcB / log(10)
    lcB <- -2.5 * log(lcB, base = 10)
  }

  mu <- theta[1]
  sigma <- theta[2]
  tau <- theta[3]

  # sorting time given delta
  time.d <- time - delta
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  time.comb <- time.temp[ord]
  leng.time.comb <- length(time.comb)
  
  # indicator taking on 1 for X(t - delta) and 0 for X(t)
  ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]

  if (micro == 0) {
    mat.temp <- matrix(c(rep(1, leng.time)), ncol = 1)
  } else if (micro == 1) {
    mat.temp <- matrix(c(rep(1, leng.time), time.d), ncol = 2)
  } else if (micro == 2) {
    mat.temp <- matrix(c(rep(1, leng.time), time.d, time.d^2), ncol = 3)
  } else if (micro == 3) {
    mat.temp <- matrix(c(rep(1, leng.time), time.d, time.d^2, time.d^3), ncol = 4)
  }

  c.pred <- mat.temp %*% c

  lc.temp <- c(lcA, lcB - c.pred)
  lc.comb <- lc.temp[ord]
  se.lc.temp <- c(se.lcA, se.lcB)
  se.lc.comb <- se.lc.temp[ord]
   
  # a.i, i = 1, ..., 2n
  time.comb.diff <- diff(time.comb)
  a.i <- exp( -time.comb.diff / tau)

  # x.i, i = 1, 2, ..., 2n
  X <- X - mu
  x <- lc.comb - mu

  # a.i, i = 2, ..., 2n
  time.comb.diff <- diff(time.comb)
  a.i <- exp( -time.comb.diff / tau )
  
  # B, i = 1, 2, ..., 2n to be saved
  B <- rep(NA, leng.time.comb)

  # mu.i, i = 1, 2, ..., 2n to be saved
  mu.i <- rep(NA, leng.time.comb)
  
  # shrinkages
  var0 <- tau * sigma^2 / 2
  B[1] <- se.lc.comb[1]^2 / (se.lc.comb[1]^2 + var0 * (1 - a.i[1]^2))
  mu.i[1] <- (1 - B[1]) * x[1] + B[1] * a.i[1] * X[2]
  X[1] <- rnorm(1, mean = mu.i[1], sd = sqrt(se.lc.comb[1]^2 * (1 - B[1])))

  for (k in 2 : (leng.time.comb - 1)) {
    B[k] <- se.lc.comb[k]^2 / 
              ( se.lc.comb[k]^2 + var0 * (1 - a.i[k - 1]^2) * (1 - a.i[k]^2) / (1 - (a.i[k - 1] * a.i[k])^2) )
    mu.i[k] <- (1 - B[k]) * x[k] + 
               B[k] * (a.i[k] * (1 - a.i[k - 1]^2) * X[k + 1] + a.i[k - 1] * (1 - a.i[k]^2) * X[k - 1]) / 
                      (1 - (a.i[k - 1] * a.i[k])^2)
    X[k] <- rnorm(1, mean = mu.i[k], sd = sqrt(se.lc.comb[k]^2 * (1 - B[k])))
  }

  B[leng.time.comb] <- se.lc.comb[leng.time.comb]^2 / 
                       (se.lc.comb[leng.time.comb]^2 + var0 * (1 - a.i[leng.time.comb - 1]^2))
  mu.i[leng.time.comb] <- (1 - B[leng.time.comb]) * x[leng.time.comb] + 
                          B[leng.time.comb] * a.i[leng.time.comb - 1] * X[leng.time.comb - 1]
  X[leng.time.comb] <- rnorm(1, mean = mu.i[leng.time.comb], 
                                sd = sqrt(se.lc.comb[leng.time.comb]^2 * (1 - B[leng.time.comb])))

  X + mu

}     



### Bayesian approach to time delay estimation
bayesian <- function(data, data.flux, 
                     theta.ini, 
                     delta.ini, delta.uniform.range, delta.proposal.scale, 
                     tau.proposal.scale, tau.prior.shape, tau.prior.scale, 
                     sigma.prior.shape, sigma.prior.scale,                        
                     asis = TRUE, micro, 
                     adaptive.freqeuncy = 100,
                     adaptive.delta = TRUE, adaptive.delta.factor = 0.01,
                     adaptive.tau = TRUE, adaptive.tau.factor = 0.01,
                     sample.size = 50, warmingup.size = 50) {

  if (asis == TRUE & adaptive.delta == TRUE & adaptive.tau == TRUE) {
    print("Options for the Bayesian method: ASIS for beta, Adaptive MCMC for Delta, and Adaptive MCMC for tau")
  } else if (asis == TRUE & adaptive.delta == TRUE & adaptive.tau == FALSE) {
    print("Options for the Bayesian method: ASIS for beta and Adaptive MCMC for Delta")
  } else if (asis == TRUE & adaptive.delta == FALSE & adaptive.tau == TRUE) {
    print("Options for the Bayesian method: ASIS for beta and Adaptive MCMC for tau")
  } else if (asis == FALSE & adaptive.delta == TRUE & adaptive.tau == TRUE) {
    print("Options for the Bayesian method:  Adaptive MCMC for Delta and tau")
  } else if (asis == TRUE & adaptive.delta == FALSE & adaptive.tau == FALSE) {
    print("Options for the Bayesian method: ASIS for beta")
  } else if (asis == FALSE & adaptive.delta == TRUE & adaptive.tau == FALSE) {
    print("Options for the Bayesian method:  Adaptive MCMC for Delta")
  } else if (asis == FALSE & adaptive.delta == FALSE & adaptive.tau == TRUE) {
    print("Options for the Bayesian method:  Adaptive MCMC for tau")
  } else if (asis == FALSE & adaptive.delta == FALSE & adaptive.tau == FALSE) {
    print("Options for the Bayesian method: None")
  }

  print(paste("Starting time:", Sys.time()))

  total.sample.size <- sample.size + warmingup.size

  time <- data[, 1]
  leng.time <- length(time)

  lcA <- data[, 2]
  se.lcA <- data[, 3]
  lcB <- data[, 4]
  se.lcB <- data[, 5]

  if (data.flux == TRUE) {
    # transform into log scale
    se.lcA <- se.lcA * 2.5 / lcA / log(10)
    lcA <- -2.5 * log(lcA, base = 10)
    se.lcB <- se.lcB * 2.5 / lcB / log(10)
    lcB <- -2.5 * log(lcB, base = 10)
  }

  delta.t <- delta.ini  # delta ini
  mu.t <- theta.ini[1]  # mu ini
  sigma.t <- theta.ini[2]  #sigma ini
  tau.t <- theta.ini[3]  #tau ini

  # sorting time given delta
  time.d <- time - delta.t
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  ti2 <- time.d^2
  ti3 <- time.d^3

  if (micro == 0) {
    lm.res <- lm(lcB - mean(lcA) ~ 1)
  } else if (micro == 1) {
    lm.res <- lm(lcB - mean(lcA) ~ time.d)
  } else if (micro == 2) {
    lm.res <- lm(lcB - mean(lcA) ~ time.d + ti2)
  } else if (micro == 3) {
    lm.res <- lm(lcB - mean(lcA) ~ time.d + ti2 + ti3)
  }

  resid <- lm.res$residuals
  X.t <- c(lcA, resid + mean(lcA))[ord]  # X ini
  c.t <- lm.res$coefficients  # beta ini

  mu.out <- rep(NA, total.sample.size)
  sigma.out <- rep(NA, total.sample.size)
  tau.out <- rep(NA, total.sample.size)
  delta.out <- rep(NA, total.sample.size)
  c.out <- matrix(NA, nrow = total.sample.size, ncol = (micro + 1))

  tau.accept <- rep(0, total.sample.size)
  delta.accept <- rep(0, total.sample.size)
  
  tau.jumps <- tau.proposal.scale * rnorm(total.sample.size)
  tau.thresh <- -rexp(total.sample.size)  
  tau.proposal.scale.adapt <- 1

  delta.jumps <- delta.proposal.scale * rnorm(total.sample.size)
  delta.thresh <- -rexp(total.sample.size)
  delta.proposal.scale.adapt <- 1


  for (i in 1 : total.sample.size) {

    # delta and X(t) update
    delta.p <- delta.t + delta.proposal.scale.adapt * delta.jumps[i]
    l.metrop <- logpostDelta(delta.p, data, c(mu.t, sigma.t, tau.t), c.t, 
                             log = data.flux, unif = delta.uniform.range, micro) -
	            logpostDelta(delta.t, data, c(mu.t, sigma.t, tau.t), c.t, 
                             log = data.flux, unif = delta.uniform.range, micro)

    if (l.metrop > delta.thresh[i]) { 
        delta.t <- delta.p 
        delta.accept[i] <- 1
        X.t <- postX(data, X.t, theta = c(mu.t, sigma.t, tau.t), 
                      delta = delta.t, c = c.t, log = data.flux, micro)
    }
	
    delta.out[i] <- delta.t  

    time.d <- time - delta.t
    time.temp <- c(time, time.d)
    ord <- order(time.temp)
    ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]
    leng.X <- length(X.t)
    lc.comb <- c(lcA, lcB)[ord]
    se.lc.comb <- c(se.lcA, se.lcB)[ord]
    time.sort <- time.temp[ord]

    # c update
    if (micro == 0) {
      T.mat <- matrix(c(rep(1, leng.time)), ncol = 1)
    } else if (micro == 1) {
      T.mat <- matrix(c(rep(1, leng.time), time.d), ncol = 2)
    } else if (micro == 2) {
      T.mat <- matrix(c(rep(1, leng.time), time.d, time.d^2), ncol = 3)
    } else if (micro == 3) {
      T.mat <- matrix(c(rep(1, leng.time), time.d, time.d^2, time.d^3), ncol = 4)
    }

    V.inv.mat <- diag(1 / se.lcB^2)
    c.t.var.inv.temp <- t(T.mat) %*% V.inv.mat %*% T.mat
    c.t.var.temp <- chol2inv(chol(c.t.var.inv.temp))
    c.t.mean.temp <- c.t.var.temp %*% t(T.mat) %*% V.inv.mat %*% (lcB - X.t[ind == 1])
    c.t.var <- chol2inv(chol(c.t.var.inv.temp + diag(micro + 1) / 10^5))
    c.t.mean <- c.t.var %*% c.t.var.inv.temp %*% c.t.mean.temp
    c.t <- c.out[i, ] <- rmnorm(1, mean = c.t.mean, varcov = c.t.var)


    if (asis == TRUE) {

      if (micro == 0) {
        T.mat <- matrix(c(rep(1, 2 * leng.time)), ncol = 1) 
      } else if (micro == 1) {
        T.mat <- matrix(c(rep(1, 2 * leng.time), time.sort), ncol = 2) 
      } else if (micro == 2) {
        T.mat <- matrix(c(rep(1, 2 * leng.time), time.sort, time.sort^2), ncol = 3) 
      } else if (micro == 3) {
        T.mat <- matrix(c(rep(1, 2 * leng.time), time.sort, time.sort^2, time.sort^3), ncol = 4) 
      }

      K.t <- X.t + T.mat %*% c.t * ind
      K.t.cent <- K.t - mu.t
      time.diff <- diff(time.sort)
      # i = 2, 3, ..., 2n
      a.i <- exp( - time.diff / tau.t )

      y.c <- c(K.t.cent[1], K.t.cent[-1] - a.i * K.t.cent[-leng.X])
      if (micro == 0) {
        X.c <- c(ind[1] * T.mat[1, ], ind[-1] * T.mat[-1, ] - a.i * ind[-leng.X] * T.mat[-leng.X, ])
      } else {
        X.c <- rbind(ind[1] * T.mat[1, ], 
                     ind[-1] * T.mat[-1, ] - a.i * ind[-leng.X] * T.mat[-leng.X, ])
      }

      V.inv.elem <- c(1, 1 / (1 - a.i^2)) / (tau.t * sigma.t^2 / 2)
      XVX <- t(X.c * V.inv.elem) %*% X.c
      c.var <- chol2inv(chol(t(X.c * V.inv.elem) %*% X.c + diag(micro + 1) / 10^5))
      c.mean <- c.var %*% t(X.c * V.inv.elem) %*% y.c
      c.t <- c.out[i, ] <- rmnorm(1, mean = c.mean, varcov = c.var)
      X.t <- K.t - T.mat %*% c.t * ind    # synchronization

    }

    # theta update
    tau.jump.adapt <- tau.proposal.scale.adapt * tau.jumps[i]
    theta.update <- postTheta(data, X = X.t, delta = delta.t, 
                               previous.theta = c(mu.t, sigma.t, tau.t), 
                               tau.jump = tau.jump.adapt, tau.thresh = tau.thresh[i],
                               tau.prior.a = tau.prior.shape, tau.prior.b = tau.prior.scale, 
                               sigma.prior.a =  sigma.prior.shape, sigma.prior.b = sigma.prior.scale)

    if (theta.update[3] != tau.t) {
      tau.accept[i] <- 1
    }
 
    mu.t <- mu.out[i] <- theta.update[1]
    sigma.t <- sigma.out[i] <- theta.update[2]
    tau.t <- tau.out[i] <- theta.update[3]

    if (adaptive.delta == TRUE) {
      if (i %% adaptive.freqeuncy == 0) {
        if(mean(delta.accept[i - (adaptive.freqeuncy - 1) : i]) > 0.44) {
          scale.adj <- exp(min(adaptive.delta.factor, 1 / sqrt(i / adaptive.freqeuncy)))
        } else if (mean(delta.accept[i - (adaptive.freqeuncy - 1) : i]) < 0.23) {
          scale.adj <- exp(-min(adaptive.delta.factor, 1 / sqrt(i / adaptive.freqeuncy)))
        } else {
          scale.adj <- 1
        }
        delta.proposal.scale.adapt <- delta.proposal.scale.adapt * scale.adj
      }
    }

    if (adaptive.tau == TRUE) {
      if (i %% adaptive.freqeuncy == 0) {
        if(mean(tau.accept[i - (adaptive.freqeuncy - 1) : i]) > 0.44) {
          scale.adj <- exp(min(adaptive.tau.factor, 1 / sqrt(i / adaptive.freqeuncy)))
        } else if (mean(tau.accept[i - (adaptive.freqeuncy - 1) : i]) < 0.23) {
          scale.adj <- exp(-min(adaptive.tau.factor, 1 / sqrt(i / adaptive.freqeuncy)))
        } else {
          scale.adj <- 1
        }
        tau.proposal.scale.adapt <- tau.proposal.scale.adapt * scale.adj
      }
    }

  }

  print(paste("Ending time:", Sys.time()))

  out <- list(delta = delta.out[-c(1 : warmingup.size)],
              beta = c.out[-c(1 : warmingup.size)],
              mu = mu.out[-c(1 : warmingup.size)], 
              sigma = sigma.out[-c(1 : warmingup.size)], 
              tau = tau.out[-c(1 : warmingup.size)],
              tau.accept.rate = mean(tau.accept),
              delta.accept.rate = mean(delta.accept))

  out

}

