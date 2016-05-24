`bisoreg` <- function(x, y, m=floor(n/2), priors=list(), inits=list(),
                      n.sim=5000, n.burn=1000, n.thin=10,
                      seed){

  ## Set up data
  ## Standardize the x's
  n <- length(y)
  m <- floor(m)
  xmin <- min(x)
  xrng <- diff(range(x))
  xstd <- (x - xmin)/xrng
  ## Set up the W matrix
  W <- matrix(unlist(lapply(1:m, function(k) pbeta(xstd, k, m - k + 1))), ncol=m)

  if (!missing(seed)) set.seed(seed)

  ## Initial values
  invals <- list(sigsq=rgamma(1,1,1),
                 tausq=rgamma(1,1,1),
                 p=rbeta(1,1,1),
                 u0=rnorm(1,0,10),
                 gam=rbinom(m,1,0.5))
  invals$u <- invals$gam*abs(rnorm(m,0,1))
  invals[names(inits)] <- inits

  ## Prior values
  privals <- list(a.sig=0.1, b.sig=0.1,
                  a.tau=1, b.tau=1,
                  a.p=1, b.p=1,
                  m0=0, ssq0=1000)
  privals[names(priors)] <- priors

  mcmctime <- system.time(gibbsout <- .C("isogibbs",
                                         ## MCMC parameters
                                         as.integer(n.sim*n.thin + n.burn),
                                         as.integer(n.thin),
                                         as.integer(n.burn),
                                         ## Constants/Data
                                         as.integer(n),
                                         as.integer(m),
                                         as.double(y),
                                         as.double(as.numeric(t(W))),
                                         ## Parameters: u[1], ..., u[n]
                                         u=as.double(invals$u),
                                         ## Parameters: tausq
                                         tausq=as.double(invals$tausq),
                                         as.double(privals$a.tau),
                                         as.double(privals$b.tau),
                                         ## Parameter: p
                                         p=as.double(invals$p),
                                         as.double(privals$a.p),
                                         as.double(privals$b.p),
                                         ## Parameter: sigsq
                                         sigsq=as.double(invals$sigsq),
                                         a.sig=as.double(privals$a.sig),
                                         b.sig=as.double(privals$b.sig),
                                         ## Parameter: u0
                                         u0=as.double(invals$u0),
                                         as.double(privals$m0),
                                         as.double(privals$ssq0),
                                         ## Parameter: gamma
                                         gam=as.integer(invals$gam),
                                         ## Results storage
                                         udraws=double(n.sim*m),
                                         tausqdraws=double(n.sim),
                                         pdraws=double(n.sim),
                                         sigsqdraws=double(n.sim),
                                         u0draws=double(n.sim),
                                         gamdraws=integer(n.sim*m),
                                         dev=double(n.sim),
                                         PACKAGE="bisoreg"))

  cat("Elapsed time: ", round(mcmctime[3]), " seconds.\n", sep="")

  u     <- matrix(gibbsout$udraws, ncol=m, byrow=TRUE)
  u0    <- gibbsout$u0draws
  tausq <- gibbsout$tausqdraws
  sigsq <- gibbsout$sigsqdraws
  gam   <- matrix(gibbsout$gamdraws, ncol=m, byrow=TRUE)
  p     <- gibbsout$pdraws
  dev   <- gibbsout$dev

  ## Compute DIC ##
  Dbar  <- mean(dev)
  mn    <- mean(u0) + W%*%apply(u, 2, mean)
  stdev <- sqrt(mean(sigsq))
  Dhat  <- -2*sum(dnorm(y, mn, stdev, log=TRUE))
  pD    <- Dbar - Dhat
  pD2   <- 0.5*var(dev)
  DIC   <- Dbar + pD
  DIC2  <- Dbar + pD2

  out       <- list()
  out$mcmctime <- mcmctime
  out$postdraws <- data.frame(u=u, u0=u0, tausq=tausq, sigsq=sigsq, gam=gam, p=p, dev=dev)
  out$W     <- W
  out$x     <- x
  out$y     <- y
  out$xstd  <- xstd
  out$xmin  <- xmin
  out$xrng  <- xrng
  out$m     <- m
  out$DIC   <- DIC
  out$DIC2  <- DIC2
  out$pD    <- pD
  out$pD2   <- pD2
  attr(out,"class") <- "biso"
  return(out)
  }


as.mcmc.biso <- function(x, ...){
  as.mcmc(x$postdraws)
  }


summary.biso <- function(object, ...){
  nms <- names(object$postdraws)
  idraws <- object$postdraws[grep("^gam", nms)]
  ddraws <- object$postdraws[grep("^[^\\(gam\\)]", nms)]
  prob0 <- c(colMeans(idraws), rep(NA, 5))
  s <- t(rbind(colMeans(ddraws), sapply(ddraws, fivenum), prob0))
  colnames(s) <- c("Mean", "Min", "1stQ", "Median", "3rdQ", "Max", "Prob(u)=0")
  return(s)
  }


print.biso <- function(x, digits=3, ...){
  s <- summary(x)
  print(s, digits=digits, ...)
  }
