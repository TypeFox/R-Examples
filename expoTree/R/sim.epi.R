sim.epi <- function(N,beta,mu,psi,max.samples,
                    min.outbreak=min(10,max.samples),max.time=-1.0) 
{
  pars <- c(N,beta,mu,psi)
  out <- .Call("RSimEpi",parameters=pars,max_samples=max.samples,
               min_outbreak=min.outbreak,max_time=max.time)
  o <- order(out$times)
  out$times <- out$times[o]
  out$ttypes <- out$ttypes[o]
  out
}


