pois = function(obs, sim, ...) {
  nlogLike = -sum(obs*log(sim) - sim, na.rm=TRUE)
  return(nlogLike)
}

normp = function(obs, sim, ...) {
  penalty = sum((sim)^2, na.rm=TRUE)
  return(penalty)
}

norm2 = function(obs, sim, ...) {
  nlogLike = sum((obs-sim)^2, na.rm=TRUE)
  return(nlogLike)
}

lnorm2 = function(obs, sim, tiny=1e-2, ...) {
  if(all(!is.finite(sim))) return(Inf)
  obs = log(obs + tiny)
  sim = log(sim + tiny)
  nlogLike = sum((obs-sim)^2, na.rm=TRUE)
  return(nlogLike)
}

lnorm3  = function(obs, sim, tiny = 1e-2, ...) {
  if(all(!is.finite(sim))) return(Inf)
  ratio = obs/sim
  ratio[is.nan(ratio)] = NA
  q = mean(ratio, na.rm=TRUE)
  obs = log(obs+tiny) 
  sim = log(sim+tiny)
  nlogLike = sum((obs-sim-log(q))^2, na.rm=TRUE)
  return(nlogLike)
}


multinom = function(sim, obs, size=20) {
  A = ncol(sim)
  sim.sum = rowSums(sim, na.rm=TRUE)
  obs.sum = rowSums(obs, na.rm=TRUE)
  sim.sum[sim.sum==0] = NA
  obs.sum[obs.sum==0] = NA  
  
  Psim     = sim/sim.sum
  Pobs     = obs/obs.sum
  
  sigma2 = ((1-Pobs)*Pobs + 1/A)/size
  
  error = log(exp(-((Pobs - Psim)^2)/(2*sigma2) + 0.001))
  
  nlogLike = -size*sum(error, na.rm=TRUE)
  return(nlogLike)
}
