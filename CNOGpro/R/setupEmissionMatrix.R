setupEmissionMatrix <-
function(nstates, mean, variance, absmax, includeZeroState=T,errorRate){
  
  emission <- matrix(NA, nrow=nstates, ncol=(absmax+1))
  
  for (state in 1:nstates){
    # Converting to Poisson-Gamma mixture distribution parameterized through p and r
    pr <- solvePR(mean*state, variance*state)
    p <- pr[1]
    r <- pr[2]
    # Setting the probability according to the probability mass function of the Poisson-Gamma compound distribution:
    for (obs in 0:absmax){
      emission[state,(obs+1)] <- calculateProb(p,r,obs)
    }
  }
  
  rownames(emission) <- seq(1,nstates)
  colnames(emission) <- seq(0,absmax)
  
  if(includeZeroState){
    # The probability mass function for state=0 is modeled using the geometric distribution with p=(1-epsilon), where epsilon is the
    # estimated fraction of erroneously mapped reads.
    zerorow <- dgeom(x=0:absmax, prob=(1-errorRate))
    emission <- rbind(zerorow, emission)
    rownames(emission) <- seq(0,nstates)
  }
  
  return(emission)
}
