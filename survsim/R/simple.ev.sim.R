simple.ev.sim <-
function(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z=NULL, beta=0, eff=0, 
           dist.ev, dist.cens, i)
  {
    nid          <- vector()
    start        <- vector()  
    stop         <- vector()
    obs          <- vector()
    it           <- vector()
    tb           <- vector()
    az1          <- NA
    time         <- vector()
    a.ev         <- NA
    b.ev         <- NA
    a.cens       <- NA
    b.cens       <- NA
    
    obs[1]   <- 1
    k.ev     <- 1
    sum      <- 0
    if (!is.null(z) && z[[1]][1] == "gamma")    az1 <- rgamma(1, as.numeric(z[[1]][2]), as.numeric(z[[1]][3]))
    if (!is.null(z) && z[[1]][1] == "exp")      az1 <- rgamma(1, 1, as.numeric(z[[1]][2]))
    if (!is.null(z) && z[[1]][1] == "weibull")  az1 <- rweibull(1, as.numeric(z[[1]][2]), as.numeric(z[[1]][3]))
    if (!is.null(z) && z[[1]][1] == "unif")     az1 <- runif(1, as.numeric(z[[1]][2]), as.numeric(z[[1]][3]))
    if (!is.null(z) && z[[1]][1] == "invgauss") az1 <- rinvgauss(1, as.numeric(z[[1]][2]), as.numeric(z[[1]][3]))
    if (is.null(z))                             az1 <- 1

    # Time to censorship
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
    start[1]   <- 0
    k.ev       <- 1
    nid[1]     <- i
    if (k.ev > length(beta0.ev))     k.ev   <- length(beta0.ev)
    
    suma <- 0
    if (!is.na(beta[1])) suma <- sum(sapply(beta, "[", k.ev) * eff)
    if (dist.ev[k.ev] == 'llogistic')
    {
      tb[1] <- az1*exp(rlogis(1, beta0.ev[k.ev] + suma, anc.ev[k.ev]))
    }else{
      if (dist.ev[k.ev] == 'weibull')
      {
        a.ev   <- anc.ev[k.ev]
        b.ev   <- (1/exp(-anc.ev[k.ev]*(beta0.ev[k.ev] + suma)))^(1/anc.ev[k.ev])
        tb[1]  <- az1*rweibull(1, a.ev, b.ev)
      }else{
        if (dist.ev[k.ev] == 'lnorm')
        {
          tb[1]  <- az1*rlnorm(1, beta0.ev[k.ev] + suma, anc.ev[k.ev])
        } #if
      } #if
    } #if  
    it[1]      <- 0
    time[1]    <- tc
    if (tb[1] < tc)
    {
      it[1]   <- 1
      time[1] <- tb[1]
    }
    
    stop[1]  <-  time[1]
      
    if (start[1] < foltime && stop[1] > foltime)
    {
      stop[1]   <- foltime
      time[1]   <- foltime
      it[1]     <- 0
    }
    
    if (start[1] < 0 && stop[1] > 0)
    {
      start[1]  <- 0
      time[1]   <- stop[1]
    }
    
    sim.ind <- data.frame(nid=nid, status=it, start=start, stop=stop, z=az1)
    for (i in 1:length(eff))
    {
      sim.ind <- cbind(sim.ind, x = eff[i])
    }
    sim.ind <- subset(sim.ind, start < foltime & stop > 0)
    
    return(sim.ind)
  }
