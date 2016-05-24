regmixEM.lambda = function (y, x, lambda = NULL, beta = NULL, sigma = NULL,
  k = 2, addintercept = TRUE, arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, 
  maxit = 10000, verb = FALSE) {
  if (arbmean == FALSE && arbvar == FALSE) {
    stop(paste("Must change constraints on beta and/or sigma!", 
               "\n"))
  }
  s = sigma
  if (addintercept) {
    x = cbind(1, x)
  }
  n <- length(y)
  p <- ncol(x)
  tmp <- regmix.lambda.init(y = y, x = x, lambda = lambda, beta = beta, 
                            s = s, k = k, addintercept = addintercept,
                            arbmean = arbmean, arbvar = arbvar)
  lambda <- tmp$lambda
  beta <- tmp$beta
  s <- tmp$s
  k <- tmp$k
  diff <- 1
  iter <- 0
  xbeta <- x %*% beta
  res2 <- (y - xbeta)^2
  if (arbmean == FALSE) {
    res2 <- sapply(1:k, function(i) res2)
  }
  comp <- t(t(lambda)/sqrt(2 * pi * s^2)) * exp(-t(t(res2)/(2 * s^2)))
  obsloglik <- sum(log(apply(comp, 1, sum)))
  ll <- obsloglik
  while (diff > epsilon && iter < 1) {
    V=as.double(sweep(lambda, 2, s+rep(0,k), "/"))
    W=as.double(sweep(res2, 2, 2*(s+rep(0,k))^2, "/"))
    z <- matrix(.C("newz", as.integer(n), as.integer(k), V=V, W=W,
                   newz=double(n*k), PACKAGE = "mixtools")$newz, ncol=k)
	z = z/apply(z,1,sum)    
    if (addintercept) {
      lm.out <- lapply(1:k, function(i) lm(y ~ x[,-1], weights = z[, i]))
    }
    else lm.out <- lapply(1:k, function(i) lm(y ~ x - 
                                              1, weights = z[, i]))
    beta.new <- sapply(lm.out, coef)
    if (arbmean == FALSE) {
      beta.new <- apply(t(apply(z, 2, sum) * t(beta.new)), 
                        1, sum)/n
    }
    xbeta.new <- x %*% beta.new
    res2 <- (y - xbeta.new)^2
    if (arbmean == FALSE) {
      res2 <- sapply(1:k, function(i) res2)
    }
    if (arbvar) {
      s.new <- sqrt(sapply(1:k, function(i)
                           sum(z[,i] * (res2[, i]))/sum(z[, i])))
    }
    else s.new <- sqrt(sum(z * res2)/n)
    beta <- beta.new
    xbeta <- x %*% beta
    s <- s.new
    sing <- sum(s < 1e-08)
    comp <- lapply(1:k, function(i)
                   lambda[,i] * dnorm(y, xbeta[, i * arbmean +
                                               (1 - arbmean)],
                                      s[i * arbvar + 
                                        (1 - arbvar)]))
    comp <- sapply(comp, cbind)
    compsum <- apply(comp, 1, sum)
    newobsloglik <- sum(log(compsum))
    if (newobsloglik < obsloglik || sing > 0 || abs(newobsloglik) == 
        Inf || is.nan(newobsloglik)){# || sum(z) != n) {
      cat("Need new starting values due to singularity...", 
          "\n")
      tmp <- regmix.lambda.init(y = y, x = x, k = k,
                                addintercept = addintercept, 
                                arbmean = arbmean, arbvar = arbvar)
      lambda <- tmp$lambda
      beta <- tmp$beta
      s <- tmp$s
      k <- tmp$k
      diff <- 1
      iter <- 0
      xbeta <- x %*% beta
      res2 <- (y - xbeta)^2
      if (arbmean == FALSE) {
        res2 <- sapply(1:k, function(i) res2)
      }
      comp <- t(t(lambda)/sqrt(2 * pi * s^2)) * exp(-t(t(res2)/(2 * 
                                                                s^2)))
      obsloglik <- sum(log(apply(comp, 1, sum)))
      ll <- obsloglik
    }
    else {
      diff <- newobsloglik - obsloglik
      obsloglik <- newobsloglik
      ll <- c(ll, obsloglik)
      iter <- iter + 1
    }
  }
  scale.order = order(s)
  sigma.min = min(s)
  if (arbmean == FALSE) {
    z = z[, scale.order]
    names(beta) <- c(paste("beta", ".", 0:(p - 1), sep = ""))
    colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a = list(x = x, y = y, lambda = lambda[,scale.order], 
      beta = beta, sigma = sigma.min, scale = s[scale.order]/sigma.min, 
      loglik = obsloglik, posterior = z[, scale.order], 
      all.loglik = ll, ft="regmixEM.lambda")
    class(a) = "mixEM"
    a
  }
  else {
    rownames(beta) <- c(paste("beta", ".", 0:(p - 1), sep = ""))
    colnames(beta) <- c(paste("comp", ".", 1:k, sep = ""))
    colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a = list(x = x, y = y, lambda = lambda, beta = beta, 
      sigma = s, loglik = obsloglik, posterior = z, all.loglik = ll,
      ft="regmixEM.lambda")
    class(a) = "mixEM"
    a
  }
}
