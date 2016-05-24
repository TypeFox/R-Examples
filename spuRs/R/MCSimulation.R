# program spuRs/resources/scripts/MCSimulation.r
# loadable spuRs function

# Discrete time Markov Chain Simulation
# This function simulates a discrete MC with transition matrix P,
# state space {0,1,..,n}, and initial state i for nsteps transitions

MCSimulation<-function(P, i, nsteps){
  n <- nrow(P)-1
  statehist <- c(i, rep(0, nsteps))
  for (step in 2:(nsteps+1)){
    statehist[step] <- sample(0:n,1,prob=P[1+statehist[step-1],])
  }
  return(statehist)
}
