
########################################################
########  THE MAIN MCMC FUNCTION FOR ABBA METHOD #######
########################################################

abba <- function(y, sites, iters, Qb = NULL, knots = sites, X = cbind(1,sites), 
  beta = NULL, alpha = 0.5, logbw = 0, tau = c(1,1,1),
  logs = matrix(0, nrow=nf, ncol=nt), u = matrix(0.5, nrow=nf, ncol=nt),
  MHbeta = matrix(rep(c(0.15,0.03,0.015),each = n), ncol=3), MHalpha = 0.01, MHlogbw = 0, 
  MHs = matrix(0.5, nrow=nf, ncol=nt), MHu = matrix(2.5, nrow=nf, ncol=nt),
  pribeta = c(10,10,10), prialpha = c(1,1), prilogbw = c(0,1), pritau = c(0.1,0.1,0.1),
  trace = 0) {
  
  # LOAD PACKAGE
  
  # BULLET PROOFING
  
  if(!is.matrix(y) || !is.numeric(y)) 
    stop("y should be a numeric matrix")
  n <- nrow(y)
  nt <- ncol(y)
  if(!is.matrix(sites) || ncol(sites) != 2) 
    stop("sites should be a matrix with two columns")
  if(nrow(sites) != n) 
    stop("number of rows in y and sites do not match")
  if(!is.matrix(knots) || ncol(knots) != 2) 
    stop("knots should be a matrix with two columns")
  nf <- nrow(knots)
  
  if(!is.list(X)) {
    if(!is.matrix(X) || nrow(X) != n)
      stop("number of rows in y and X do not match")
    X <- list(mu = X, lsig = X, xi = X)
  } else {
    if(length(X) != 3)
      stop("X should be a list of length three")
    if(nrow(X[[1]]) != n)
      stop("number of rows in y and X[[1]] do not match")
    if(nrow(X[[2]]) != n)
      stop("number of rows in y and X[[2]] do not match")
    if(nrow(X[[3]]) != n)
      stop("number of rows in y and X[[3]] do not match")
  }
  p <- sapply(X, ncol)
  
  if(!is.null(Qb)) {
    if(!is.matrix(Qb) || nrow(Qb) != ncol(Qb))
      stop("Qb must be a square matrix")
    if(nrow(Qb) != n)
      stop("number of rows in y and Qb do not match")
    if(!isSymmetric(Qb))
      stop("Qb must be symmetric")
  } else {
    Qb <- -1/fields::rdist(sites,sites)^2
    diag(Qb) <- 0
    diag(Qb) <- -rowSums(Qb)
    Qb <- -Qb/min(Qb)
  }
  
  if(length(alpha) != 1)
    stop("alpha must be a single value")
  if(alpha <= 0 || alpha >= 1)
    stop("alpha must be in the interval (0,1)")
  if(length(logbw) != 1)
    stop("logbw must be a single value")
  if(length(tau) != 3)
    stop("tau must be of length three")
  if(any(tau <= 0))
    stop("tau values must be positive")
  if(!is.matrix(logs) || !is.matrix(u))
    stop("logs and u must be matrices")
  if(nrow(logs) != nf || nrow(u) != nf)
    stop("number of rows of knots, logs and u must all match")
  if(ncol(logs) != nt || ncol(u) != nt)
    stop("number of columns of y, logs and u must all match")
  if(any(u <= 0) || any(u >= 1))
    stop("u values must be in the interval (0,1)")
  
  if(!is.matrix(MHbeta) || !is.matrix(MHs) || !is.matrix(MHu))
    stop("MHbeta, MHs and MHu must be matrices")
  if(nrow(MHbeta) != n) 
    stop("number of rows in y and MHbeta do not match")
  if(ncol(MHbeta) != 3) 
    stop("MHbeta should have three columns")
  if(nrow(MHs) != nf || nrow(MHu) != nf)
    stop("number of rows of knots, MHs and MHu must all match")
  if(ncol(MHs) != nt || ncol(MHu) != nt)
    stop("number of columns of y, MHs and MHu must all match")
  if(any(MHbeta < 0) || any(MHs < 0) || any(MHu < 0))
    stop("MHbeta, MHs and MHu values must be non-negative")
  if(length(MHalpha) != 1 || length(MHlogbw) != 1)
    stop("MHalpha and MHlogbw must be single values")
  if(MHalpha < 0 || MHlogbw < 0)
    stop("MHalpha and MHlogbw must be non-negative")
  
  if(length(pribeta) != 3 || any(pribeta <= 0))
    stop("pribeta must be a positive vector of length three")
  if(length(prialpha) != 2 || any(prialpha <= 0))
    stop("prialpha must be a positive vector of length two")
  if(length(prilogbw) != 2 || prilogbw[2] <= 0)
    stop("prilogbw must be a vector of length two")
  if(length(pritau) != 3 || any(pritau <= 0))
    stop("pritau must be a positive vector of length three")
  if(length(trace) != 1 || trace < 0)
    stop("trace must be a positive integer")
  
  if(any(is.na(sites)))
    stop("site locations cannot have missing values")
  if(any(is.na(knots)))
    stop("knot locations cannot have missing values")
  if(any(apply(y, 1, function(x) all(is.na(x)))))
    stop("a site cannot have all missing values")

  dw2 <- fields::rdist(sites,knots)^2
  dw2[dw2 < 0.0001] <- 0
  
  # INITIAL VALUES
  
  if(!is.null(beta)) {
    if(!is.matrix(beta) || ncol(beta) != 3)
      stop("beta must be a matrix with three columns")
    if(nrow(beta) != n)
      stop("number of rows in y and beta do not match")
  } else {
    beta <- matrix(0, nrow=n, ncol=3)
    for(j in 1:n) {
      # pwm estimates assuming zero shape
      yj <- y[j,][!is.na(y[j,])]
      if(length(yj) == 1) 
        stop("cannot compute starting values with one data point")
      lambda0 <- mean((1:length(yj) - 1)/(length(yj) - 1) * sort(yj))
      sig0 <- (2*lambda0 - mean(yj))/log(2)
      mu0 <- mean(yj) + digamma(1) * sig0
      beta[j,] <- c(mu0, log(sig0), 0)
    }
  }
  mnb <- list(mu = rep(0, p[1]), lsig = rep(0, p[2]), xi = rep(0, p[3]))
  for(l in 1:3) mnb[[l]][1] <- colMeans(beta)[l]
  
  # C INTERFACE PREPARATION 
  
  iters <- as.integer(iters)
  ns <- as.integer(n)
  nf <- as.integer(nf)
  nt <- as.integer(nt)
  
  y <- as.double(t(y))
  npmu <- ncol(X[[1]])
  Xmu <- as.double(t(X[[1]]))
  npsig <- ncol(X[[2]])
  Xsig <- as.double(t(X[[2]]))
  npxi <- ncol(X[[3]])
  Xxi <- as.double(t(X[[3]]))
  
  mnbmu <- as.double(mnb[[1]])
  taumu <- as.double(tau[1])
  mnbsig <- as.double(mnb[[2]])
  tausig <- as.double(tau[2])
  mnbxi <- as.double(mnb[[3]])
  tauxi <- as.double(tau[3])
  
  alpha <- as.double(alpha)
  logbw <- as.double(logbw)
  
  beta <- as.double(t(beta))
  logs <- as.double(t(logs))
  u <- as.double(t(u))
  
  Qb <- as.double(t(Qb))
  dw2 <- as.double(t(dw2))
  
  MHbeta <- as.double(t(MHbeta))
  MHalpha <- as.double(MHalpha)
  MHlogbw <- as.double(MHlogbw)
  MHs <- as.double(t(MHs))
  MHu <- as.double(t(MHu))
  
  pribeta <- as.double(pribeta)
  prialpha <- as.double(prialpha)
  prilogbw <- as.double(prilogbw)
  pritau <- as.double(pritau)
  
  trace <- as.integer(trace)
  
  # C INTERFACE 
  
  np <- npmu + npsig + npxi + 6
  mc <- .C("fgevspatial", iters, ns, nf, nt, y, Xmu, npmu, Xsig, npsig, Xxi, npxi,
    mnbmu, taumu, mnbsig, tausig, mnbxi, tauxi,
    alpha, logbw, beta, logs, u, Qb, dw2,
    MHbeta, MHalpha, MHlogbw, MHs, MHu,
    pribeta, prialpha, prilogbw, pritau, trace,     
    betasamp = double(iters * ns * 3),  
    paramsamp = double(iters * np),
    psrvsamp = double(iters * nf * nt),
    urvsamp = double(iters * nf * nt), NAOK = TRUE)
  
  # RETURN OUTPUT
  mc$paramsamp <- matrix(mc$paramsamp, nrow=iters, ncol=np, byrow=TRUE)
  nms <- paste(rep(c("loc","lscale","shape"), p+1), c(paste0("b",0:(p[1]-1)),"tau", 
    paste0("b",0:(p[2]-1)),"tau", paste0("b",0:(p[3]-1)),"tau"), sep=".")
  rownames(mc$paramsamp) <- paste0("i", 1:iters)
  colnames(mc$paramsamp) <- c(nms,"alpha","logbw","logpost")
  
  mc$betasamp <- aperm(array(mc$betasamp, dim = c(3,ns,iters)), 3:1)
  dimnames(mc$betasamp) <- list(paste0("i", 1:iters), paste0("s", 1:n), c("loc","lscale","shape"))
  mc$psrvsamp <- aperm(array(mc$psrvsamp, dim = c(nt,nf,iters)), 3:1)
  dimnames(mc$psrvsamp) <- list(paste0("i", 1:iters), paste0("k", 1:nf), paste0("t", 1:nt))
  mc$urvsamp <- aperm(array(mc$urvsamp, dim = c(nt,nf,iters)), 3:1)
  dimnames(mc$urvsamp) <- list(paste0("i", 1:iters), paste0("k", 1:nf), paste0("t", 1:nt))
  
  list(beta.samples = mc$betasamp, param.samples = mc$paramsamp, 
       psrv.samples = mc$psrvsamp, urv.samples = mc$urvsamp)
}

##########################################################
########  THE MAIN MCMC FUNCTION FOR LATENT METHOD #######
##########################################################

abba_latent <- function(y, sites, iters, Qb = NULL, X = cbind(1,sites), 
  beta = NULL, tau = c(1,1,1),
  MHbeta = matrix(rep(c(0.15,0.03,0.015),each = n), ncol=3), 
  pribeta = c(10,10,10), pritau = c(0.1,0.1,0.1),
  trace = 0) {
  
  # BULLET PROOFING
  
  if(!is.matrix(y) || !is.numeric(y)) 
    stop("y should be a numeric matrix")
  n <- nrow(y)
  nt <- ncol(y)
  if(!is.matrix(sites) || ncol(sites) != 2) 
    stop("sites should be a matrix with two columns")
  if(nrow(sites) != n) 
    stop("number of rows in y and sites do not match")
  
  if(!is.list(X)) {
    if(!is.matrix(X) || nrow(X) != n)
      stop("number of rows in y and X do not match")
    X <- list(mu = X, lsig = X, xi = X)
  } else {
    if(length(X) != 3)
      stop("X should be a list of length three")
    if(nrow(X[[1]]) != n)
      stop("number of rows in y and X[[1]] do not match")
    if(nrow(X[[2]]) != n)
      stop("number of rows in y and X[[2]] do not match")
    if(nrow(X[[3]]) != n)
      stop("number of rows in y and X[[3]] do not match")
  }
  p <- sapply(X, ncol)
  
  if(!is.null(Qb)) {
    if(!is.matrix(Qb) || nrow(Qb) != ncol(Qb))
      stop("Qb must be a square matrix")
    if(nrow(Qb) != n)
      stop("number of rows in y and Qb do not match")
    if(!isSymmetric(Qb))
      stop("Qb must be symmetric")
  } else {
    Qb <- -1/as.matrix(dist(sites))^2
    diag(Qb) <- 0
    diag(Qb) <- -rowSums(Qb)
    Qb <- -Qb/min(Qb)
  }
  
  if(length(tau) != 3)
    stop("tau must be of length three")
  if(any(tau <= 0))
    stop("tau values must be positive")
  
  if(!is.matrix(MHbeta))
    stop("MHbeta must be a matrix")
  if(nrow(MHbeta) != n) 
    stop("number of rows in y and MHbeta do not match")
  if(ncol(MHbeta) != 3) 
    stop("MHbeta should have three columns")
  if(any(MHbeta < 0))
    stop("MHbeta values must be non-negative")
  
  if(length(pribeta) != 3 || any(pribeta <= 0))
    stop("pribeta must be a positive vector of length three")
  if(length(pritau) != 3 || any(pritau <= 0))
    stop("pritau must be a positive vector of length three")
  if(length(trace) != 1 || trace < 0)
    stop("trace must be a positive integer")
  
  if(any(is.na(sites)))
    stop("site locations cannot have missing values")
  if(any(apply(y, 1, function(x) all(is.na(x)))))
    stop("a site cannot have all missing values")
  
  # INITIAL VALUES
  
  if(!is.null(beta)) {
    if(!is.matrix(beta) || ncol(beta) != 3)
      stop("beta must be a matrix with three columns")
    if(nrow(beta) != n)
      stop("number of rows in y and beta do not match")
  } else {
    beta <- matrix(0, nrow=n, ncol=3)
    for(j in 1:n) {
      # pwm estimates assuming zero shape
      yj <- y[j,][!is.na(y[j,])]
      if(length(yj) == 1) 
        stop("cannot compute starting values with one data point")
      lambda0 <- mean((1:length(yj) - 1)/(length(yj) - 1) * sort(yj))
      sig0 <- (2*lambda0 - mean(yj))/log(2)
      mu0 <- mean(yj) + digamma(1) * sig0
      beta[j,] <- c(mu0, log(sig0), 0) 
    }
  }
  mnb <- list(mu = rep(0, p[1]), lsig = rep(0, p[2]), xi = rep(0, p[3]))
  for(l in 1:3) mnb[[l]][1] <- colMeans(beta)[l]
  
  # C INTERFACE PREPARATION 
  
  iters <- as.integer(iters)
  ns <- as.integer(n)
  nt <- as.integer(nt)
  
  y <- as.double(t(y))
  npmu <- ncol(X[[1]])
  Xmu <- as.double(t(X[[1]]))
  npsig <- ncol(X[[2]])
  Xsig <- as.double(t(X[[2]]))
  npxi <- ncol(X[[3]])
  Xxi <- as.double(t(X[[3]]))
  
  mnbmu <- as.double(mnb[[1]])
  taumu <- as.double(tau[1])
  mnbsig <- as.double(mnb[[2]])
  tausig <- as.double(tau[2])
  mnbxi <- as.double(mnb[[3]])
  tauxi <- as.double(tau[3])
  
  beta <- as.double(t(beta))
  Qb <- as.double(t(Qb))
  MHbeta <- as.double(t(MHbeta))
  pribeta <- as.double(pribeta)
  pritau <- as.double(pritau)
  trace <- as.integer(trace)
  
  # C INTERFACE 
  
  np <- npmu + npsig + npxi + 4
  mc <- .C("fgevspatial_latent", iters, ns, nt, y, Xmu, npmu, Xsig, npsig, Xxi, npxi,
           mnbmu, taumu, mnbsig, tausig, mnbxi, tauxi,
           beta, Qb, MHbeta, 
           pribeta, pritau, trace,     
           betasamp = double(iters * ns * 3),  
           paramsamp = double(iters * np), NAOK = TRUE)
  
  # RETURN OUTPUT
  mc$paramsamp <- matrix(mc$paramsamp, nrow=iters, ncol=np, byrow=TRUE)
  nms <- paste(rep(c("loc","lscale","shape"), p+1), c(paste0("b",0:(p[1]-1)),"tau", 
                                                      paste0("b",0:(p[2]-1)),"tau", paste0("b",0:(p[3]-1)),"tau"), sep=".")
  rownames(mc$paramsamp) <- paste0("i", 1:iters)
  colnames(mc$paramsamp) <- c(nms,"logpost")
  
  mc$betasamp <- aperm(array(mc$betasamp, dim = c(3,ns,iters)), 3:1)
  dimnames(mc$betasamp) <- list(paste0("i", 1:iters), paste0("s", 1:n), c("loc","lscale","shape"))
 
  list(beta.samples = mc$betasamp, param.samples = mc$paramsamp)
}
