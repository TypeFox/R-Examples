BiCopGofKendall <- function(u1, u2, family, B = 100, level = 0.05) {
  
  if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
    stop("u1 and/or u2 are not set or have length zero.")
  if (any(u1 > 1) || any(u1 < 0)) 
    stop("Data has be in the interval [0,1].")
  if (any(u2 > 1) || any(u2 < 0)) 
    stop("Data has be in the interval [0,1].")
  if (length(u1) != length(u2)) 
    stop("Lengths of 'u1' and 'u2' do not match.")
  if (length(u1) < 2) 
    stop("Number of observations has to be at least 2.")
  if (!(family %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20, 23, 24, 26, 27, 28, 
    29, 30, 33, 34, 36, 37, 38, 39, 40))) 
    stop("Copula family not implemented.")
  if (level < 0 & level > 1) 
    stop("Significance level has to be between 0 and 1.")
  
  if (family %in% c(13, 14, 16, 17, 18, 19, 20)) {
    u1 <- 1 - u1
    u2 <- 1 - u2
    family <- family - 10
  } else if (family %in% c(23, 24, 26, 27, 28, 29, 30)) {
    u1 <- 1 - u1
    family <- family - 20
  } else if (family %in% c(33, 34, 36, 37, 38, 39, 40)) {
    u2 <- 1 - u2
    family <- family - 30
  }
  
  ostat <- obs.stat(u1, u2, family)
  
  if (B == 0) {
    # no bootstrap
    
    sn.obs <- ostat$Sn
    tn.obs <- ostat$Tn
    out <- list(Sn = sn.obs, Tn = tn.obs)
    
  } else {
    bstat <- list()
    for (i in 1:B) bstat[[i]] <- boot.stat(u1, u2, family)
    
    sn.boot <- rep(0, B)
    tn.boot <- rep(0, B)
    for (i in 1:B) {
      sn.boot[i] <- bstat[[i]]$sn
      tn.boot[i] <- bstat[[i]]$tn
    }
    
    sn.obs <- ostat$Sn
    tn.obs <- ostat$Tn
    
    k <- as.integer((1 - level) * B)
    sn.critical <- sn.boot[k]  # critical value of test at level 0.05
    tn.critical <- tn.boot[k]  # critical value of test at level 0.05
    
    
    pv.sn <- sapply(sn.obs, function(x) (1/B) * length(which(sn.boot[1:B] >= x)))  # P-value of Sn
    pv.tn <- sapply(tn.obs, function(x) (1/B) * length(which(tn.boot[1:B] >= x)))  # P-value of Tn
    
    out <- list(p.value.CvM = pv.sn, p.value.KS = pv.tn, statistic.CvM = sn.obs, statistic.KS = tn.obs)
    
  }
  
  return(out)
}


boot.stat <- function(u, v, fam) {
  
  n <- length(u)
  t <- seq(1, n)/(n + 1e-04)
  kt <- rep(0, n)
  
  # estimate paramemter for different copula family from (u,v)
  param <- suppressWarnings({
    BiCopEst(u, v, family = fam)
  })
  # calulate k(t) and kn(t) of bootstrap sample data
  if (fam == 1) {
    # normal
    sam <- BiCopSim(n, 1, param$par, param$par2)  # generate data for the simulation of K(t)
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # parameter estimation of sample data
    sim <- BiCopSim(10000, 1, sam.par$par, sam.par$par2)  # generate data for the simulation of theo. K(t)
    cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
    dcop <- rep(0, 10000)
    for (i in 1:10000) dcop[i] <- pmvnorm(upper = c(qnorm(sim[i, 1]), qnorm(sim[i, 2])), corr = cormat)
    kt <- sapply(t, function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
  }
  if (fam == 2) {
    # t
    sam <- BiCopSim(n, fam, param$par, param$par2)  # generate data for the simulation of K(t)
    sam.par <- suppressWarnings({
      BiCopEst(sam[, 1], sam[, 2], family = fam)
    })  # parameter estimation of sample data
    sim <- BiCopSim(10000, fam, sam.par$par, sam.par$par2)  # generate data for the simulation of theo. K(t)
    cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
    dcop <- rep(0, 10000)
    for (i in 1:10000) dcop[i] <- pmvt(upper = c(qt(sim[i, 1], df = param$par2), qt(sim[i, 2], df = param$par2)), 
      corr = cormat, df = param$par2)
    kt <- sapply(t, function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
  } else if (fam == 3) {
    # Clayton
    sam <- BiCopSim(n, 3, param$par)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data
    kt <- t + t * (1 - t^sam.par)/sam.par
  } else if (fam == 4) {
    # gumbel
    sam <- BiCopSim(n, 4, param$par)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data
    kt <- t - t * log(t)/(sam.par)
  } else if (fam == 5) {
    # frank
    sam <- BiCopSim(n, 5, param$par)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data
    kt <- t + log((1 - exp(-sam.par))/(1 - exp(-sam.par * t))) * (1 - exp(-sam.par * t))/(sam.par * exp(-sam.par * 
      t))
  } else if (fam == 6) {
    # joe
    sam <- BiCopSim(n, 6, param$par)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)$par  # estimate parameter of sample data        
    kt <- t - (log(1 - (1 - t)^sam.par) * (1 - (1 - t))^sam.par)/(sam.par * (1 - t)^(sam.par - 1))
  } else if (fam == 7) {
    # BB1
    sam <- BiCopSim(n, 7, param$par, param$par2)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
    theta <- sam.par$par
    delta <- sam.par$par2
    kt <- t + 1/(theta * delta) * (t^(-theta) - 1)/(t^(-1 - theta))
  } else if (fam == 8) {
    # BB6
    sam <- BiCopSim(n, 8, param$par, param$par2)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
    theta <- sam.par$par
    delta <- sam.par$par2
    kt <- t + log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * theta)
  } else if (fam == 9) {
    # BB7
    sam <- BiCopSim(n, 9, param$par, param$par2)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
    theta <- sam.par$par
    delta <- sam.par$par2
    kt <- t + 1/(theta * delta) * ((1 - (1 - t)^theta)^(-delta) - 1)/((1 - t)^(theta - 1) * (1 - (1 - 
      t)^theta)^(-delta - 1))
  } else if (fam == 10) {
    # BB8
    sam <- BiCopSim(n, 10, param$par, param$par2)  # generate sample data
    sam.par <- BiCopEst(sam[, 1], sam[, 2], family = fam)  # estimate parameter of sample data
    theta <- sam.par$par
    delta <- sam.par$par2
    kt <- t + log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * delta)^(-theta) + 
      (1 - t * delta)^(-theta) * t * delta)/(theta * delta)
  }
  
  # calculate emp. Kn
  w <- rep(0, n)
  w[1:n] <- mapply(function(x, y) (1/n) * length(which(x > sam[, 1] & y > sam[, 2])), sam[, 1], sam[, 2])
  w <- sort(w)
  kn <- rep(0, n)
  kn <- sapply(t, function(x) (1/n) * length(which(w[1:n] <= x)))
  
  # calculate test statistic Sn
  Sn1 <- 0
  Sn2 <- 0
  for (j in 1:(n - 1)) {
    Sn1 <- Sn1 + ((kn[j])^2 * (kt[j + 1] - kt[j]))
    Sn2 <- Sn2 + (kn[j]) * ((kt[j + 1])^2 - (kt[j])^2)
  }
  sn <- n/3 + n * Sn1 - n * Sn2
  # calculation of test statistics Tn
  tm <- matrix(0, n - 1, 2)
  # mit i=0
  for (j in 1:(n - 1)) {
    tm[j, 1] <- abs(kn[j] - kt[j])
  }
  # mit i=1
  for (j in 1:(n - 1)) {
    tm[j, 2] <- abs(kn[j] - kt[j + 1])
  }
  tn <- max(tm) * sqrt(n)
  
  sn <- sort(sn)  # vector of ordered statistic Sn
  tn <- sort(tn)  # vector of ordered statistic Tn
  out <- list(sn = sn, tn = tn)
}



obs.stat <- function(u, v, fam) {
  
  n <- length(u)
  t <- seq(1, n)/(n + 1e-04)
  kt <- rep(0, n)
  
  # estimate paramemter for different copula family from (u,v)
  param <- suppressWarnings({
    BiCopEst(u, v, family = fam)
  })
  
  # calculate observed K(t) of (u,v)
  kt.obs <- rep(0, n)
  if (fam == 1) {
    sim <- BiCopSim(10000, 1, param$par)  # generate data for the simulation of K(t)
    cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
    dcop <- rep(0, 10000)
    for (i in 1:10000) dcop[i] <- pmvnorm(upper = c(qnorm(sim[i, 1]), qnorm(sim[i, 2])), corr = cormat)
    kt.obs <- sapply(t, function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
  } else if (fam == 2) {
    sim <- BiCopSim(10000, 2, param$par, param$par2)  # generate data for the simulation of K(t)
    cormat <- matrix(c(1, param$par, param$par, 1), 2, 2)
    dcop <- rep(0, 10000)
    for (i in 1:10000) dcop[i] <- pmvt(upper = c(qt(sim[i, 1], df = param$par2), qt(sim[i, 2], df = param$par2)), 
      corr = cormat, df = param$par2)
    kt.obs <- sapply(t, function(x) (1/10000) * length(which(dcop[1:10000] <= x)))  # simulate K(t) of sample data
  } else if (fam == 3) {
    kt.obs <- t + t * (1 - t^param$par)/param$par
  } else if (fam == 4) {
    kt.obs <- t - t * log(t)/(param$par)
  } else if (fam == 5) {
    kt.obs <- t + log((1 - exp(-param$par))/(1 - exp(-param$par * t))) * (1 - exp(-param$par * t))/(param$par * 
      exp(-param$par * t))
  } else if (fam == 6) {
    kt.obs <- t - (log(1 - (1 - t)^param$par) * (1 - (1 - t))^param$par)/(param$par * (1 - t)^(param$par - 
      1))
  } else if (fam == 7) {
    theta <- param$par
    delta <- param$par2
    kt.obs <- t + 1/(theta * delta) * (t^(-theta) - 1)/(t^(-1 - theta))
  } else if (fam == 8) {
    theta <- param$par
    delta <- param$par2
    kt.obs <- t + log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * 
      theta)
  } else if (fam == 9) {
    theta <- param$par
    delta <- param$par2
    kt.obs <- t + 1/(theta * delta) * ((1 - (1 - t)^theta)^(-delta) - 1)/((1 - t)^(theta - 1) * (1 - 
      (1 - t)^theta)^(-delta - 1))
  } else if (fam == 10) {
    theta <- param$par
    delta <- param$par2
    kt.obs <- t + log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * 
      delta)^(-theta) + (1 - t * delta)^(-theta) * t * delta)/(theta * delta)
  }
  # calculation of observed Kn
  w <- rep(0, n)
  w[1:n] <- mapply(function(x, y) (1/n) * length(which(x > u & y > v)), u, v)
  
  w <- sort(w)
  kn.obs <- rep(0, n)
  kn.obs <- sapply(t, function(x) (1/n) * length(which(w[1:n] <= x)))
  
  # calculation of observed value Sn
  Sn1 <- 0
  Sn2 <- 0
  for (j in 1:(n - 1)) {
    Sn1 <- Sn1 + ((kn.obs[j])^2 * (kt.obs[j + 1] - kt.obs[j]))
    Sn2 <- Sn2 + (kn.obs[j]) * ((kt.obs[j + 1])^2 - (kt.obs[j])^2)
  }
  Sn <- n/3 + n * Sn1 - n * Sn2  # observed Sn
  
  # calculation of observed Tn
  tn.obs <- matrix(0, n - 1, 2)
  # mit i=0
  for (j in 1:(n - 1)) {
    tn.obs[j, 1] <- abs(kn.obs[j] - kt.obs[j])
  }
  # mit i=1
  for (j in 1:(n - 1)) {
    tn.obs[j, 2] <- abs(kn.obs[j] - kt.obs[j + 1])
  }
  Tn <- max(tn.obs) * sqrt(n)
  out <- list(Sn = Sn, Tn = Tn)
  return(out)
} 
