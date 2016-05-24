"epi.indirectadj" <- function(obs, pop, std, units, conf.level = 0.95){   
    # How many strata (rows) are there?
    n.strata <- dim(pop)[1]
        
    # How many covariates are there?
    n.cov <- dim(pop)[2]

    N <- 1 - ((1 - conf.level) / 2)
    alpha <- 1 - conf.level
    z <- qnorm(N, mean = 0, sd = 1)
    
    tmp <- data.frame(strata = rep(rownames(pop), times = n.cov), cov = rep(colnames(pop), each = n.strata), 
       pop = as.vector(pop), std = rep(as.vector(std[1:n.cov]), each = n.strata))

    # Expected events (equals std incidence multiplied by population size):
    tmp$exp <- (tmp$pop * tmp$std)
    # tmp <- tmp[order(tmp$strata, tmp$cov),] 
    
    # Crude risk by strata:
    # Turn 'obs' into a table object so calculations can easily be done by strata:
    t.obs <- by(data = obs, INDICES = rownames(obs), FUN = sum)
    t.exp <- by(data = tmp$exp, INDICES = tmp$strata, FUN = sum)
    t.pop <- by(data = tmp$pop, INDICES = tmp$strata, FUN = sum)
    
    # Confidence interval for crude incidence risk estimates corrected following email from Gillian Raab:
    crude.p <- t.obs / t.pop
    # crude.se <- crude.p / sqrt(t.pop)                            ## Incorrect.
    crude.se <- crude.p / sqrt(t.obs)                              ## replaced pop by obs
    crude.l <- qchisq(alpha / 2, 2 * t.obs) / 2 / t.pop            ## next 2 lines changed
    crude.u <- qchisq(1 - alpha / 2, 2 *(t.obs + 1)) / 2 / t.pop
    crude.strata <- data.frame(est = as.vector(crude.p) * units, lower = as.vector(crude.l) * units, 
       upper = as.vector(crude.u) * units)
    rownames(crude.strata) <- names(t.exp)
    
    # Indirectly adjusted risk for each strata (see page 378 of Stata manual):
    t.obs <- by(data = obs, INDICES = rownames(obs), FUN = sum)
    t.exp <- by(data = tmp$exp, INDICES = tmp$strata, FUN = sum)
    t.pop <- by(data = tmp$pop, INDICES = tmp$strata, FUN = sum)
    
    if(n.cov > 1){
       adj.p <- (std[n.cov + 1] * (t.obs / t.exp))
       adj.l <- (std[n.cov + 1] * (qpois((1 - conf.level) / 2, lambda = t.obs, log.p = FALSE) / t.exp))
       adj.u <- (std[n.cov + 1] * (qpois(1 - (1 - conf.level) / 2, lambda = t.obs, log.p = FALSE) / t.exp))
       adj.strata <- data.frame(est = as.vector(adj.p) * units, lower = as.vector(adj.l) * units, upper = as.vector(adj.u) * units)
       rownames(adj.strata) <- names(t.exp)
       }

   if(n.cov == 1){
      adj.p <- (std * (t.obs / t.exp))
      adj.l <- (std * (qpois((1 - conf.level) / 2, lambda = t.obs, log.p = FALSE) / t.exp))
      adj.u <- (std * (qpois(1 - (1 - conf.level) / 2, lambda = t.obs, log.p = FALSE) / t.exp))
      adj.strata <- data.frame(est = as.vector(adj.p) * units, lower = as.vector(adj.l) * units, upper = as.vector(adj.u) * units)
      rownames(adj.strata) <- names(t.exp)
      }

    # Crude standardised mortality ratio (confidence intervals based on Breslow and Day 1987 p 69-71): 
    smr.p <- t.obs / t.exp
    smr.l <- qpois((1 - conf.level) / 2, lambda = t.obs, log.p = FALSE) / t.exp
    smr.u <- qpois(1 - (1 - conf.level) / 2, lambda = t.obs, log.p = FALSE) / t.exp
    smr.strata <- data.frame(obs = as.vector(t.obs), exp = as.vector(t.exp), est = as.vector(smr.p), lower = as.vector(smr.l), upper = as.vector(smr.u))
    rownames(smr.strata) <- names(t.exp)   
    
    rval <- list(crude.strata = crude.strata, adj.strata = adj.strata, smr.strata = smr.strata)
    return(rval)  
}