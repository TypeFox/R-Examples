setupTransitionMatrix <-
function(nstates, remainprob, includeZeroState=T){
  
  if(includeZeroState){
    nstates <- nstates+1
  }
  changeprob <- 1 - remainprob
  per_state_prob <- changeprob/(nstates-1)
  transition <- matrix(per_state_prob, nrow=nstates, ncol=nstates)
  
  for (i in 1:nstates){
    transition[i,i] <- remainprob
  }
  
  return(transition)
}
