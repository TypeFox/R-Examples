###################################################
### chunk number 1: 
###################################################

# Programme to simulate epidemies which were
# introduced by point sources.
# The basis of this proagramme is a combination of
# a Hidden Markov Modell (to get random dates
# for outbreaks) and a simple Model to simulate
# the epidemy.
#
# Parameters:
# r - probability to get a new epidemy at time i if there was one
#     at time i-1
# p - probability to get no new epidemy at time i if there was none
#     at time i-1
# length - number of timesteps to visit
#
# Parameters for the background:
# A - Amplitude, default = 1.
# alpha - Incidence, default = 1.
# beta - time dependent regression coefficient, default = 0.
# phi - weeks of seaonal move, default = 0.
# frequency - frequency of the sinus, default = 1.
# state - a eventually given markov chain,
#               which defines the status at this time (outbreak or not)
# K - additional weigth for an outbreak

sim.pointSource <- function(p = 0.99, r = 0.01, length = 400, A = 1, alpha = 1, beta = 0,
                                phi = 0, frequency = 1, state = NULL, K){

  if(is.null(state)){
        # create a markov-chain
    state <- matrix(data = 0, ncol = 1, nrow = length)
    state[1] <- 0 #hoehle - fix: rbinom(1,1,0.5) # always begin with a zero

        # create the transition matrix
    transitionMatrix <- matrix(data = c(p, (1-r),(1-p), r), nrow = 2, ncol = 2)

    if(length(state) > 1){ # just do it if there is a preceding value
      for (i in 2:length){
                        # check the matrix for the correct line and take the right
                        # probability. The last value of state is the newest.
        state[i] <- rbinom(1,1,transitionMatrix[state[i-1] + 1, 2])
      }
    }
  }

  # go sure to have the rigth length as parameter
  length <- length(state)
  observed <-sim.seasonalNoise(A, alpha, beta, phi, length, frequency, state, K)$seasonalBackground

  result <- list(observed = observed, state = state, A = A, alpha = alpha, beta = beta, K = K, p = p, r = r, freq=52, start=c(2001,1))
  class(result) = "disProg" # for disease progress

  return(result)
}





