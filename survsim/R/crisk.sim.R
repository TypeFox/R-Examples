crisk.sim <-
  function (n, foltime, dist.ev, anc.ev, beta0.ev, dist.cens = "weibull", 
            anc.cens, beta0.cens, z = NULL, beta = NA, x = NA, nsit) 
  {
    # Arguments check
    if (length(anc.ev) != length(beta0.ev)) stop("Wrong number of parameters")
    if (length(anc.cens) != length(beta0.cens) && length(anc.cens) != 1) stop("Wrong number of parameters")
    if (length(anc.ev) != length(dist.ev)) stop("Wrong number of parameters")
    if (!is.na(x) && is.na(beta)) stop("Wrong specification of covariables")
    if (is.na(x) && !is.na(beta)) stop("Wrong specification of covariables")
    if (length(beta0.ev) != nsit) stop("Wrong number of distributions or events")
    if (!is.null(z) && length(z) != nsit && length(z) != 1) stop("Wrong numbers of elements in z")
    if (!is.null(z) && !all(lapply(z, function(x) x[1]) %in% c("unif","weibull","invgauss", "gamma","exp"))) 
      stop("Wrong specification of z")
    if (!is.null(z) && any(lapply(z, function(x) length(x)) != 3))
    {
      if(any(lapply(z[lapply(z, function(x) length(x)) != 3], function(x) length(x)) != 2)) stop("Wrong specification of z")
      if(any(lapply(z[lapply(z, function(x) length(x)) != 3], function(x) length(x)) == 2))
      {
        for (i in 1:length(z[lapply(z, function(x) length(x)) == 2]))
        {
          if (z[lapply(z, function(x) length(x)) != 3][[i]][1] != "exp") stop("Wrong specification of z")
        }
      }
    }
    if(!is.null(z) && any(lapply(z, function(x) x[1]) == "unif"))
    {
      for (i in 1:length(z[lapply(z, function(x) x[1]) == "unif"]))
      {
        if (as.numeric(z[lapply(z, function(x) x[1]) == "unif"][[i]][2])-as.numeric(z[lapply(z, function(x) x[1]) == "unif"][[i]][3]) >= 0) 
           stop("Wrong specification of z")
        if (as.numeric(z[lapply(z, function(x) x[1]) == "unif"][[i]][2]) < 0) stop("Wrong specification of z")
      }
    }
    
    sim.data <- list()
    eff <- vector()
    eff[1] <- 0
    for (i in 1:n) {
      if (!is.na(x[1])) {
        for (k in 1:length(x)) {
          if (x[[k]][1] == "unif") 
            eff[k] <- runif(1, as.numeric(x[[k]][2]), as.numeric(x[[k]][3]))
          if (x[[k]][1] == "normal") 
            eff[k] <- rnorm(1, as.numeric(x[[k]][2]), as.numeric(x[[k]][3]))
          if (x[[k]][1] == "bern") 
            eff[k] <- rbinom(1, 1, as.numeric(x[[k]][2]))
        }
      }
      sim.data[[i]] <- crisk.ncens.sim(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z, beta, eff, 
                                       dist.ev, dist.cens, i, nsit)
    }
    sim.data <- do.call(rbind, sim.data)
    sim.data$cause[sim.data$status==0] <- NA
    class(sim.data) <- c("crisk.data.sim", "data.frame")
    attr(sim.data, "n") <- n
    attr(sim.data, "foltime") <- foltime
    attr(sim.data, "nsit") <- nsit
    return(sim.data)
  }