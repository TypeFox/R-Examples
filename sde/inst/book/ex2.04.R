# ex2.04.R
BM.1 <- function(N=10000){ # brutal code
 W <- NULL 
 for(i in 2:(N+1))
        W <- c(W, rnorm(1) / sqrt(N))	
}

BM.2 <- function(N=10000){ # smarter
 W <- numeric(N+1) 
 Z <- rnorm(N)
 for(i in 2:(N+1))
        W[i] <- W[i-1] + Z[i-1] / sqrt(N)	
}

BM.vec <- function(N=10000) # awesome !
 W <- c(0,cumsum(rnorm(N)/sqrt(N)))

set.seed(123)
system.time(BM.1())
 
set.seed(123)
system.time(BM.2())
 
set.seed(123)
system.time(BM.vec())

