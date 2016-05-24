logit <- function(x){
  log(x/(1-x))
}

plogit <- function(x, m, s){
  pnorm(log(x/(1-x)), m ,s)
}

qlogit <- function(x, m, s){
  z <- qnorm(x, m, s)
  exp(z) / (1 + exp(z))
}

dlogit <- function(x, m, s){
  1 / (x * (1 - x)) * dnorm(log(x / (1 - x)), m, s)
}




psample <- function(medianfit, precisionfit, lower = NA, upper = NA, 
                    median.dist, precision.dist, n.rep = 10000, n.X = 100){
  
  mediandist <- getmediandist(medianfit, median.dist)
  
  f <- getdists(precisionfit$transform)
  lim <- getlimits(lower, upper, f, mediandist,  precisionfit)
  
  X <- seq(from = lim$lower, to = lim$upper, length = n.X)
  Xmat <- matrix(X, n.rep, n.X, byrow=T)
  
  mu <- matrix(mediandist$rand(n.rep, mediandist$m, mediandist$s), n.rep, n.X)
  
  if(precision.dist == "gamma"){
    sigma <- matrix(sqrt(1 / rgamma(n.rep, precisionfit$Gamma[[1]], 
                                    precisionfit$Gamma[[2]])),
                    n.rep, n.X)
  }
  
  if(precision.dist == "lognormal"){
    sigma <- matrix(sqrt(1 / rlnorm(n.rep, precisionfit$Log.normal[[1]], 
                                    precisionfit$Log.normal[[2]])),
                    n.rep, n.X)
  }
  
  pX <- f$cdf(Xmat, f$trans(mu), sigma)
  
  list(X=X, pX=pX)
}



taildensities <- function(m, s, tails, n.x, lower, upper, dens, quan, trans){
  xl <- seq(from = lower, to = quan(tails/2, m, s), 
            length = n.x)
  dl <- dens(xl, m, s)
  xu <- seq(from = quan(1-tails/2, m, s), to = upper,
            length = n.x)
  du <- dens(xu, m, s)
  data.frame(xl = xl, dl = dl, xu = xu, du = du)
}

getdists <- function(transform){
  if (transform == "identity"){
    dens <- dnorm
    quan <- qnorm
    cdf <- pnorm
    trans <- identity
  }
  
  if (transform == "log"){
    dens <- dlnorm
    quan <- qlnorm
    cdf <- plnorm
    trans <- log
  }
  
  if (transform == "logit"){
    dens <- dlogit
    quan <- qlogit
    cdf <- plogit
    trans <- logit
  }
  
  list(dens = dens, quan = quan, trans = trans, cdf = cdf)
}

getlimits <- function(lower, upper, f, mediandist, precisionfit){
  
  a<-precisionfit$Gamma[[1]]
  b<-precisionfit$Gamma[[2]]
  
  if(is.na(lower)) lower <- f$quan(0.001, 
                                   f$trans(mediandist$quan(0.001, mediandist$m, mediandist$s)),
                                   1/qgamma(0.001, a, b)^0.5)
  if(is.na(upper)) upper <- f$quan(0.999,
                                   f$trans(mediandist$quan(0.999, mediandist$m, mediandist$s)),
                                   1/qgamma(0.001, a, b)^0.5)
  list(lower = lower, upper = upper)
}

getmediandist <- function(medianfit, d){
  if(d == "best"){
    ssq <- medianfit$ssq
    ssq[is.na(ssq)] <- Inf
    if(ssq[1,1] < ssq[1,4]){d <- "normal"}else{d <- "lognormal"}
  }
  
  if(d == "normal"){
    rand <- rnorm
    quan <- qnorm
    m <- medianfit$Normal[[1]]
    s <- medianfit$Normal[[2]]
  }
  
  if(d == "lognormal"){
    rand <- rlnorm
    quan <- qlnorm
    m <- medianfit$Log.normal[[1]]
    s <- medianfit$Log.normal[[2]]
  }
  list(rand = rand, quan = quan, m = m, s=s)
}