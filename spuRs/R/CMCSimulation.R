# program spuRs/resources/scripts/CMCSimulation.r
# loadable spuRs function

# Continuous time Markov Chain Simulation
# This function simulates a continuous time MC with
# state space {0, 1, ..., n}, rate matrix Q, and initial state i
# over the period (0, Tend).  It returns the realisation
# cbind(statehist, timehist) where the vector statehist gives the 
# successive states of the jump process and the vector
# timehist the jump times.
# It also plots the realisation if plotflag is TRUE

CMCSimulation <- function(Q, i, Tend, plotflag = FALSE){
  n <- nrow(Q) - 1
  rates <- -diag(Q)
  P <- diag(1/rates) %*% Q + diag(1, n+1)
  statehist <- c(i)
  timehist <- c(0)
  time <- 0
  currentstate <- i
  jump <- 2
  
  while (time < Tend){
  	time <- time + rexp(1, rates[currentstate+1])
  	timehist[jump] <- time
  	currentstate <- sample(0:n, 1, prob=P[currentstate+1,])
  	statehist[jump] <- currentstate
  	jump <- jump + 1
  }
  
  if (plotflag) {
    plot(timehist, statehist, type="s", xlab="Time", ylab="State", 
         ylim=c(0,n), xlim=c(0,Tend), yaxt="n")
    axis(2, at=0:n)
  }
  return(cbind(statehist, timehist)) 
}