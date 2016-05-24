PTOnePeriodPD <- function(portf.uncond, portf.def, conf.interval = 0.9) {
  # One Period PD estimation Pluto & Tasche model
  # Args:
  #   portf.uncond:   unconditional portfolio distribution from the worst to the best credit quality
  #   portf.def:      number of defaults in a given rating class
  #   conf.interval:  condifence interval of PD estimates
  # Returns:
  #                   estimated conditional PDs
  r.num <- length(portf.uncond) 
  r.PD <- rep(0, r.num)
  portf.CNum <- rev(cumsum(portf.uncond))
  portf.CDef <- rev(cumsum(portf.def)) 
  
  for (r in seq_len(r.num)) {
    if (portf.CDef[r] == portf.CNum[r]) {
      r.PD[r] <- 1
    } else {
      f <- function(x) pbinom(portf.CDef[r], portf.CNum[r], x) - 1 + conf.interval
      r.PD[r] <- uniroot(f, c(0, 1))$root
    }
  }
  return(rev(r.PD))
}


PTLinkFunc <- function(pd, rho, y) {
  # Is used in multiperiod Pluto & Tasche model for conditional PD estimation (given unconditional PD,
  # correlation with systematic factor and systematic factor realization)
  # Args:
  #   pd:             unconditional PD
  #   rho:            correlation with systematic factor
  #   y:              realization of systematic factor
  # Returns:
  #                   estimated conditional PD
  return(pnorm((qnorm(pd) - sqrt(rho) * y) / (sqrt(1 - rho))))
}

PTProbLessKdef <- function(portf.uncond, portf.def, pd, rho, rSt) {
  # Estimates probability of occurrence of less than portf.def defualts given PD for Multi-period PT model
  # Args:
  #   portf.uncond:   unconditional portfolio distribution from the worst to the best credit quality
  #   portf.def:      number of defaults in a given rating class
  #   rho:            correlation with systematic factor
  #   rSt:            realizations of systematic factor (scenarios x periods)
  # Returns:
  #                   mean probability across simulations of occurrence of less than portf.def defaults given pd  
  
  piSt.onePer <- 1 - PTLinkFunc(pd, rho, rSt)
  piSt <- 1 - apply(piSt.onePer, 1, prod)
  s = 0
  for (i in seq.int(0, portf.def)) {
    s <- s + choose(portf.uncond, i) * (piSt^i) * ((1 - piSt)^(portf.uncond - i))      
  }
  return(mean(s))
}


PTMultiPeriodPD <- function(portf.uncond, portf.def, rho, cor.St, kT, kNS = 1000, conf.interval = 0.9) {
  # Estimates PDs according to multi-period Pluto & Tasche model
  # Args:
  #   portf.uncond:   unconditional portfolio distribution from the worst to the best credit quality
  #   portf.def:      number of defaults in a given rating class
  #   rho:            correlation with systematic factor
  #   cor.rSt:        correlation matrix of systematic factor realization through the time. In case constant is given - power matrix is used
  #   kT:             number of periods used in the PD estimation
  #   kNS:            number of simulations for integral estimation (using Monte-Carlo approach)
  #   conf.interval:  confidence interval
  # Returns:
  #                   conditional PDs according to multi-period Pluto & Tasche model  
  if (length(cor.St) == 1) {
    cor.ST <- matrix(numeric(kT * kT), nrow = kT) #Correlation matrix of systematic factors
    for (i in 0:(kT - 1)) 
      for (j in 1:kT)
        cor.ST[i * kT + j] <- cor.St ^ abs(i - j + 1)
  } else {
    cor.ST <- cor.St 
  }
  
  r.num <- length(portf.uncond) 
  r.PD <- rep(0, r.num)
  portf.CNum <- rev(cumsum(portf.uncond))
  portf.CDef <- rev(cumsum(portf.def)) 
  
  rSt <- MASS::mvrnorm(n = kNS, rep(0, kT), Sigma = cor.ST)
  
  for (r in seq_len(r.num)) {  # Iterating through rating classes
    if (portf.CNum[r] == portf.CDef[r]) {
      r.PD[r] <- 1
    } else {
      f <- function(x) PTProbLessKdef(portf.CNum[r], portf.CDef[r], x, rho, rSt) - 1 + conf.interval
      r.PD[r] <- uniroot(f, c(0, 1))$root
    }
  }
  return(rev(r.PD))
}
