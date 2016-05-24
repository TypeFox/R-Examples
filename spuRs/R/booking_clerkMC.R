booking_clerkMC <- function(personArrRate,
                           callArrRate,
                           personServRate,
                           callServRate,
                           t.end) {
  
# simulate the harassed booking clerk using a Markov chain
# state is a vector length 3
# first element is clerk state 0/1/2 for idle/serve person/serve call
# second element is number of people waiting
# third element is number of calls waiting

t <- 0 # start time
state <- c(0,0,0) # start state
t.hist <- t
state.hist <- state
while (t < t.end) {
  if (state[1] == 0) { #clerk idle
    rate.out <- personArrRate + callArrRate
    u <- runif(1)
    if (u < personArrRate/rate.out) { # person arrives
      state <- c(1,0,0)
    } else { # call arrives
      state <- c(2,0,0)
    }
  } else if (state[1] == 1) { # clerk serving a person
    rate.out <- personArrRate + callArrRate + personServRate
    u <- runif(1)
    if (u < personArrRate/rate.out) { # person arrives
      state[2] <- state[2] + 1
    } else if (u < (personArrRate + callArrRate)/rate.out) { 
	  # call arrives
      state[3] <- state[3] + 1
    } else { # finished serving
      if (state[2] > 0) { # still people
        state[2] <- state[2] - 1
      } else if (state[3] > 0) { # no people but calls
        state[1] <- 2
        state[3] <- state[3] - 1
      } else { # no people or calls
        state <- c(0,0,0)
      }
    }
  } else { # clerk serving a call
    rate.out <- personArrRate + callArrRate + callServRate
    u <- runif(1)
    if (u < personArrRate/rate.out) { # person arrives
      state[2] <- state[2] + 1
    } else if (u < (personArrRate + callArrRate)/rate.out) { 
	  # call arrives
      state[3] <- state[3] + 1
    } else { # finished serving
      if (state[2] > 0) { # still people
        state[1] <- 1
        state[2] <- state[2] - 1
      } else if (state[3] > 0) { # no people but calls
        state[3] <- state[3] - 1
      } else { # no people or calls
        state <- c(0,0,0)
      }
    }
  }
  t <- t + rexp(1, rate.out)
  # data recording  
  t.hist <- c(t.hist, t)
  state.hist <- rbind(state.hist, state)
}   
out <- cbind(t.hist, state.hist)
colnames(out) <- c("time", "clerk", "people", "calls")
rownames (out) <- 1:dim(out)[1]
return(out)
}