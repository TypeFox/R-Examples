booking_clerkMC <- function(person.arrival.rate,
                           call.arrival.rate,
                           person.service.rate,
                           call.service.rate,
                           t.end) {
  
# simulate the harassed booking clerk using a Markov chain
# state is a vector length 3
# first element is clerk state 0/1/2 for idle/serving a person/serving a call
# second element is number of people waiting
# third element is number of calls waiting

t <- 0 # start time
state <- c(0,0,0) # start state
t.hist <- t
state.hist <- state
while (t < t.end) {
  if (state[1] == 0) { #clerk idle
    rate.out <- person.arrival.rate + call.arrival.rate
    u <- runif(1)
    if (u < person.arrival.rate/rate.out) { # person arrives
      state <- c(1,0,0)
    } else { # call arrives
      state <- c(2,0,0)
    }
  } else if (state[1] == 1) { # clerk serving a person
    rate.out <- person.arrival.rate + call.arrival.rate + person.service.rate
    u <- runif(1)
    if (u < person.arrival.rate/rate.out) { # person arrives
      state[2] <- state[2] + 1
    } else if (u < (person.arrival.rate + call.arrival.rate)/rate.out) { # call arrives
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
    rate.out <- person.arrival.rate + call.arrival.rate + call.service.rate
    u <- runif(1)
    if (u < person.arrival.rate/rate.out) { # person arrives
      state[2] <- state[2] + 1
    } else if (u < (person.arrival.rate + call.arrival.rate)/rate.out) { # call arrives
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
return(cbind(t.hist, state.hist))
}


# sample call queue on a delta lattice from time tstart to tend
# repeat n times and average
# default rates give a heavily loaded system
call.sample <- function(tstart = 0, tend = 10000, delta = 1, n = 10, 
                        person.arrival.rate = 0.15, call.arrival.rate = 0.15) {
  m <- floor((tend-tstart)/delta) + 1
  callQ <- matrix(nrow = m, ncol = n)
  for (i in 1:n) {
    bcout <- booking_clerkMC(person.arrival.rate, call.arrival.rate, 0.3, 0.45, tend)
    callQ[,i] <- approx(bcout[,1], bcout[,4], xout=seq(tstart, tend, delta), 
                        method="constant", ties="ordered")$y
  }
  return(rowMeans(callQ))
}

# apply a moving average of width 2*w + 1
ma <- function(x, w) {
  n <- length(x)
  if (n < 2*w + 1) {
    m <- floor((n + 1)/2)
    y <- rep(0, m)
    for (i in 1:m) {
      y[i] <- mean(x[1:(2*i - 1)])
    }
  } else {
    y <- rep(0, n - w)
    for (i in 1:w) {
      y[i] <- mean(x[1:(2*i - 1)])
    }
    for (i in (w + 1):(n - w)) {
      y[i] <- mean(x[(i-w):(i+w)])
    }
  }
  return(y)
}

# effect of averaging and smoothing
# Welch's method
set.seed(20)
opar <- par(mfrow = c(3, 1))
#
tend <- 1000
delta <- 1
n <- 1
calls <- call.sample(0, tend, delta, n)
plot(seq(0, tend, delta), calls, type="l", 
     xlab="", ylab="", main=paste("n =", n))
#
n <- 10
calls <- call.sample(0, tend, delta, n)
plot(seq(0, tend, delta), calls, type="l", 
     xlab="", ylab="", main=paste("n =", n))
#     
w <- 200
callQ.w <- ma(calls, w)
plot(seq(0, tend, delta), c(callQ.w, rep(NA, w)), type="l", 
     xlab="", ylab="", main=paste("n =", n, "w =", w))
par(opar)

# choose burn in 800 and block size 400
# use blocking to estimate CI for mean queue size
set.seed(20)
burn_in <- 800
block_size <- 400 
n_blocks <- 50
delta <- 1
tend <- burn_in + block_size*n_blocks*delta
callQ <- call.sample(burn_in, tend, delta, 1)[-1]
# calculate block means
callQ_blocks <- matrix(callQ, ncol=n_blocks)
callQ_means <- colMeans(callQ_blocks)
# CI using blocks
cat("mean length call queue", mean(callQ_means), "+/-", 2*sd(callQ_means)/sqrt(n_blocks), "(using blocks)\n")
# CI without using blocks
cat("mean length call queue", mean(callQ), "+/-", 2*sd(callQ)/sqrt(length(callQ)), "(without blocks)\n")

# calculate sample autocorrelation
(rho <- acf(callQ_means, 1, plot=FALSE)$acf[2]) 
2*pnorm(abs(rho + 1/n_blocks), mean=0, sd=1/sqrt(n_blocks), lower.tail=FALSE) 

# effect of averaging and smoothing
# Welch's method
windows()
set.seed(20)
opar <- par(mfrow = c(4, 1))
#
tend <- 10000
delta <- 1
n <- 1
calls <- call.sample(0, tend, delta, n)
plot(seq(0, tend, delta), calls, type="l", 
     xlab="", ylab="", main=paste("n =", n))
#
n <- 10
calls <- call.sample(0, tend, delta, n)
plot(seq(0, tend, delta), calls, type="l", 
     xlab="", ylab="", main=paste("n =", n))
#
w <- 200
callQ.w <- ma(calls, w)
plot(seq(0, tend, delta), c(callQ.w, rep(NA, w)), type="l", 
     xlab="", ylab="", main=paste("n =", n, "w =", w))
#     
w <- 2000
callQ.w <- ma(calls, w)
plot(seq(0, tend, delta), c(callQ.w, rep(NA, w)), type="l", 
     xlab="", ylab="", main=paste("n =", n, "w =", w))
par(opar)

# choose burn in 2000 and block size 1000
# use blocking to estimate CI for mean queue size
set.seed(20)
burn_in <- 2000
block_size <- 1000 
n_blocks <- 50
delta <- 1
tend <- burn_in + block_size*n_blocks*delta
callQ <- call.sample(burn_in, tend, delta, 1)[-1]
# calculate block means
callQ_blocks <- matrix(callQ, ncol=n_blocks)
callQ_means <- colMeans(callQ_blocks)
# CI using blocks
cat("mean length call queue", mean(callQ_means), "+/-", 2*sd(callQ_means)/sqrt(n_blocks), "(using blocks)\n")
# CI without using blocks
cat("mean length call queue", mean(callQ), "+/-", 2*sd(callQ)/sqrt(length(callQ)), "(without blocks)\n")

# calculate sample autocorrelation
(rho <- acf(callQ_means, 1, plot=FALSE)$acf[2]) 
2*pnorm(abs(rho + 1/n_blocks), mean=0, sd=1/sqrt(n_blocks), lower.tail=FALSE) 


# for comparison, consider blocking for iid data
x_blocks <- matrix(rgamma(n_blocks*block_size, shape=2), ncol=n_blocks)
x_means <- apply(x_blocks, 2, mean)
cat("mean", mean(x_means), "+/-", 2*sd(x_means)/sqrt(n_blocks), "(using blocks)\n")
cat("mean", mean(x_blocks), "+/-", 2*sd(as.vector(x_blocks))/sqrt(length(x_blocks)), "(without blocks)\n")
