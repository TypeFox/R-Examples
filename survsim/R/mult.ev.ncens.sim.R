mult.ev.ncens.sim <-
function(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z=NULL, beta=0, eff=0, 
           dist.ev, dist.cens, i, nsit)
  {
    nid          <- vector()
    start        <- vector()  
    stop         <- vector()
    obs          <- vector()
    it           <- vector()
    tb           <- vector()
    az1          <- vector()
    time         <- vector()
    a.ev         <- NA
    b.ev         <- NA
    a.cens       <- NA
    b.cens       <- NA
    
    obs[1]   <- 1
    k.ev     <- 1
    sum      <- 0
    if (is.null(z))
    {
      for (l in 1:nsit)
      {
        az1[l] <- 1
      }
    }else{
      for (l in 1:length(z))
      {
        if (!is.na(z[[l]][1]) && z[[l]][1] == "gamma") 
          az1[l] <- rgamma(1, as.numeric(z[[l]][2]), as.numeric(z[[l]][3]))
        if (!is.na(z[[l]][1]) && z[[l]][1] == "exp") 
          az1[l] <- rgamma(1, 1, as.numeric(z[[l]][2]))
        if (!is.na(z[[l]][1]) && z[[l]][1] == "weibull") 
          az1[l] <- rweibull(1, as.numeric(z[[l]][2]), as.numeric(z[[l]][3]))
        if (!is.na(z[[l]][1]) && z[[l]][1] == "unif") 
          az1[l] <- runif(1, as.numeric(z[[l]][2]), as.numeric(z[[l]][3]))
        if (!is.na(z[[l]][1]) && z[[l]][1] == "invgauss") 
          az1[l] <- rinvgauss(1, as.numeric(z[[l]][2]), as.numeric(z[[l]][3]))
      }
      if (length(z) == 1)
      {
        for (l in 2:nsit)
        {
          az1[l] <- az1[1]
        }
      }
    }
    
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
    for (j in 1:nsit)
    {      
      start[j]   <- 0
      k.ev       <- j
      nid[j]     <- i
      if (k.ev > length(beta0.ev))     k.ev   <- length(beta0.ev)
      
      suma <- 0
      if (!is.na(beta[1])) suma <- sum(sapply(beta, "[", k.ev) * eff)
      if (dist.ev[k.ev] == 'llogistic')
      {
        tb[j] <- az1[j]*exp(rlogis(1, beta0.ev[k.ev] + suma, anc.ev[k.ev]))
      }else{
        if (dist.ev[k.ev] == 'weibull')
        {
          a.ev   <- anc.ev[k.ev]
          b.ev   <- (1/exp(-anc.ev[k.ev]*(beta0.ev[k.ev] + suma)))^(1/anc.ev[k.ev])
          tb[j]  <- az1[j]*rweibull(1, a.ev, b.ev)
        }else{
          if (dist.ev[k.ev] == 'lnorm')
          {
            tb[j]  <- az1[j]*rlnorm(1, beta0.ev[k.ev] + suma, anc.ev[k.ev])
          } #if
        } #if
      } #if  
      it[j]      <- 0
      time[j]    <- tc
      if (tb[j] < tc)
      {
        it[j]   <- 1
        time[j] <- tb[j]
      }
      
      stop[j]  <-  time[j]
      
      if (start[j] < foltime && stop[j] > foltime)
      {
        stop[j]   <- foltime
        time[j]   <- foltime
        it[j]     <- 0
      }
      
      if (start[j] < 0 && stop[j] > 0)
      {
        start[j]  <- 0
        time[j]   <- stop[j]
      }
      
      if (j > 1)  obs[j] <- obs[j-1]
      if (start[j] < foltime && stop[j] > 0 && j > 1) obs[j] <- obs[j-1] + 1
      
      j      <- j + 1
      k.ev   <- k.ev + 1
    } #while
    
    sim.ind <- data.frame(nid=nid, ev.num=obs, time=time, 
                          status=it, start=start, stop=stop, z=az1)
    for (i in 1:length(eff))
    {
      sim.ind <- cbind(sim.ind, x = eff[i])
    }
    sim.ind <- subset(sim.ind, start < foltime & stop > 0)
    
    return(sim.ind)
  }
