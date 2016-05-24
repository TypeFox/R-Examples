crisk.ncens.sim <- 
  function (foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z=NULL, beta=0, eff=0, 
            dist.ev, dist.cens, i, nsit) 
{
  nid    <- NA
  start  <- NA
  stop   <- NA
  obs    <- NA
  it     <- NA
  time   <- NA
  pro    <- vector()
  cause  <- NA
  a.ev   <- vector()
  b.ev   <- vector()
  a.cens <- NA
  b.cens <- NA
  obs[1] <- 1
  k.ev   <- 1
  sum    <- 0
  cshaz  <- list()
  az1    <- vector()
  
  if (is.null(z))
  {
    for (j in 1:nsit)
    {
      az1[j] <- 1
    }
  }else{
    for (j in 1:length(z))
    {
      if (!is.na(z[[j]][1]) && z[[j]][1] == "gamma") 
        az1[j] <- rgamma(1, as.numeric(z[[j]][2]), as.numeric(z[[j]][3]))
      if (!is.na(z[[j]][1]) && z[[j]][1] == "exp") 
        az1[j] <- rgamma(1, 1, as.numeric(z[[j]][2]))
      if (!is.na(z[[j]][1]) && z[[j]][1] == "weibull") 
        az1[j] <- rweibull(1, as.numeric(z[[j]][2]), as.numeric(z[[j]][3]))
      if (!is.na(z[[j]][1]) && z[[j]][1] == "unif") 
        az1[j] <- runif(1, as.numeric(z[[j]][2]), as.numeric(z[[j]][3]))
      if (!is.na(z[[j]][1]) && z[[j]][1] == "invgauss") 
        az1[j] <- rinvgauss(1, as.numeric(z[[j]][2]), as.numeric(z[[j]][3]))
    }
    if (length(z) == 1)
    {
      for (j in 2:nsit)
      {
        az1[j] <- az1[j]
      }
    }
  }
  
  if (dist.cens == "llogistic") {
    tc <- exp(rlogis(1, beta0.cens, anc.cens))
  }
  else {
    if (dist.cens == "weibull") {
      a.cens <- anc.cens
      b.cens <- (1/exp(-anc.cens * (beta0.cens)))^(1/anc.cens)
      tc <- rweibull(1, a.cens, b.cens)
    }
    else {
      if (dist.cens == "lnorm") {
        tc <- rlnorm(1, beta0.cens, anc.cens)
      }
    else {
      if (dist.cens== "unif") {
        tc <- runif(1, beta0.cens, anc.cens)
      }
    }
    }
  }
  suma <- vector()
  for (m2 in 1:nsit)
  {
    suma[m2] <- 0
    for (m1 in 1:length(beta))
    {
      suma[m2] <- suma[m2] + beta[[m1]][m2]*eff[m1]
    }
  }
  if (all(is.na(suma))) suma <- rep(0, nsit)
  # Cause-specific hazards
  for (k in 1:nsit)
  {
    if (dist.ev[k] == "llogistic") {
      a.ev[k] <- 1/exp(beta0.ev[k] + suma[k])
      b.ev[k] <- anc.ev[k]
      cshaz[[k]] <- function(t, r) {
        par1 <- eval(parse(text="a.ev[r]"))
        par2 <- eval(parse(text="b.ev[r]"))
        z    <- eval(parse(text="az1[r]"))
        return(z*(par1*par2*(t^(par2-1)))/(1+par1*(t^par2)))}
    }
    else {
      if (dist.ev[k] == "weibull") {
        a.ev[k] <- beta0.ev[k] + suma[k]
        b.ev[k] <- anc.ev[k]
        cshaz[[k]] <- function(t, r) {
          par1 <- eval(parse(text="a.ev[r]"))
          par2 <- eval(parse(text="b.ev[r]"))
          z    <- eval(parse(text="az1[r]"))
          return(z*(((1/par2)/exp(par1))^(1/par2))*t^((1/par2)-1))}
      }
      else {
        if (dist.ev[k] == "lnorm") {
          a.ev[k] <- beta0.ev[k] + suma[k]
          b.ev[k] <- anc.ev[k]
          cshaz[[k]] <- function(t, r) {
            par1 <- eval(parse(text="a.ev[r]"))
            par2 <- eval(parse(text="b.ev[r]"))
            z    <- eval(parse(text="az1[r]"))
            return(z*(dnorm((log(t)-par1)/par2)/(par2*t*(1-pnorm((log(t)-par1)/par2)))))}
        }#if
      }#if
    }#if
  }#for
  
  A <- function(t,y){ #Cumulative all-cause hazard A
    res <- 0
    for (k in 1:length(cshaz))
    {
      res <- res + integrate(cshaz[[k]], lower=0.001, upper=t, r=k, subdivisions=1000)$value
    }
    res <- res + y
    return(res[1])
  }
  u     <- runif(1)
  iters <- 0
  while (A(0.001, log(1-u))*A(foltime, log(1-u)) > 0 & iters < 1000)
  {
    u     <- runif(1)
    iters <- iters + 1
  }
  if (iters >= 1000) stop("Error: Values at endpoints not of opposite sign. \n")
  tb <- uniroot(A, c(0, foltime), tol=0.0001, y=log(1-u))$root

  sumprob <- 0
  for (k in 1:length(cshaz))
  {
    sumprob <- sumprob + cshaz[[k]](tb, k) 
  }
  for (k in 1:length(cshaz))
  {
      pro[k] <- cshaz[[k]](tb, k) / sumprob
  }
  cause1 <- rmultinom(1, 1, prob = pro)
  for (k in 1:length(cshaz))
  {
    if (cause1[k] == 1) cause <- k
  }
  az <- az1[k]
  nid <- i
  start <- 0
  it <- 0
  time <- tc
  if (tb < tc) {
      it <- 1
      time <- tb
    }
    stop <- time
    if (start < foltime && stop > foltime) {
      stop <- foltime
      time <- foltime
      it <- 0
    }
    if (start < 0 && stop > 0) {
      start <- 0
      time <- stop
    }
    
  sim.ind <- data.frame(nid = nid, cause = cause, time = time, status = it, 
                        start = start, stop = stop, z = az)
  for (k in 1:length(eff)) {
    sim.ind <- cbind(sim.ind, x = eff[k])
  }
  sim.ind <- subset(sim.ind, start < foltime & stop > 0)
  return(sim.ind)
}