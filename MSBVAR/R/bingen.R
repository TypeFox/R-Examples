# Function to generate the vector of binary indicators of the state
# space.  This is a state-space generation function for one
# observation.
#
# p = filtered probabilities of each regime for an observation
# m = number of regimes.
#
# Returns a vector of 0-1 that indicate the regime.  A one is returned
# for the regime that the observation falls into.

bingen <- function(p, Q, st1){ .Call("bingen.cpp", p, Q, st1) }

bingen.R <- function(p, Q, st1){
    h <- dim(Q)[1]
   i <- 1

    while(i<h)
    { pr0 <- p[i]*Q[st1,i]/sum(p*Q[st1,])

      if(pr0>runif(1))
      { return(diag(h)[i,]) } else { st1 <- i <- i+1 }
    }
    return(diag(h)[,h])
##     sp <- p*Q[,st1]
##     sp <- sp/sum(sp)
##     tmp <- sum(cumsum(sp) < runif(1))+ 1
##     return(diag(h)[tmp,])
}

generate.states <- function(filtered.prob, Q)
  { TT <- nrow(filtered.prob)
    h <- ncol(filtered.prob)

    # storage
    ss <- matrix(0, nrow=h, ncol=TT)

    # generate the TT state
    ss[(sum(cumsum(filtered.prob[TT,]) < runif(1)) + 1),TT] <- 1

    for (t in (TT-1):1)
      {
          ss[,t] <- bingen.R(filtered.prob[t,], Q, which(ss[,t+1]==1))
      }
    transitions <- count.transitions(t(ss))
    return(list(SS=t(ss), transitions=transitions))
  }
