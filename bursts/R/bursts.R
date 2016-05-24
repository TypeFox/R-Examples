# bursts.R

kleinberg <- function (offsets, s=2, gamma=1)
{
  if (s <= 1) {
    stop('s must be greater than 1!')
  }
  if (gamma <= 0) {
    stop('gamma must be positive!')
  }
  if (length(offsets) < 1) {
    stop('offsets must be non-empty!') 
  }
  if (length(offsets) == 1) {
    bursts <- data.frame(level=0, start=offsets[1], end=offsets[1])
    class(bursts) <- c('bursts', 'data.frame')
    return(bursts)
  }
  
  # Get the data into the right format.  The value of gaps[t] gives the length
  # of the gap between the events at offsets[t] and offsets[t+1].
  offsets <- sort(offsets)
  gaps <- offsets[2:length(offsets)] - offsets[1:length(offsets) - 1]
  gaps <- as.numeric(gaps)
  if (any(gaps == 0)) {
    stop('Input cannot contain events with zero time between!')
  }

  # Compute the average event rate, etc.
  T <- sum(gaps)
  n <- length(gaps)
  ghat <- T / n
  
  # Compute an upper bound on the number of states used in the optimal state
  # sequence, as per Kleinberg.
  k <- ceiling(1 + log(T, s) + log(1 / min(gaps), s))
  
  # Set up the transition cost function tau, and f, the probability density
  # function for gap lengths when in state j, with precomputed parameters.
  gammalogn <- gamma * log(n)
  tau <- function (i, j) if (i >= j) 0 else (j - i) * gammalogn
  alpha <- mapply(function (x) s ** x / ghat, 0:(k-1))
  f <- function (j, x) alpha[j] * exp(-alpha[j] * x)
  
  # Compute the optimal state sequence for the model, using the Viterbi
  # algorithm.  In each iteration t, we compute for each possible state j the
  # minimum costs of partial state sequences up to iteration t that end in
  # that state.  These numbers are stored in C.  We use q to keep track of
  # the optimal sequences for each j.
  C <- c(0, rep(Inf, k - 1))
  q <- matrix(NA, k, 0)
  for (t in 1:n) {
    Cprime <- rep(NA, k)
    qprime <- matrix(NA, k, t)
    for (j in 1:k) {
      cost <- mapply(function (ell) C[ell] + tau(ell, j), 1:k)
      ell <- which.min(cost)
      Cprime[j] <- cost[ell] - log(f(j, gaps[t]))
      if (t > 1)
        qprime[j, 1:(t-1)] <- q[ell, ]
      qprime[j, t] <- j
    }
    C <- Cprime
    q <- qprime
  }
  
  # Extract the state sequence with the minimum final cost.
  j <- which.min(C)
  q <- q[j,]
  
  # Compute the number of entries we will need in the output.
  prev_q <- 0
  N <- 0
  for (t in 1:n) {
    if (q[t] > prev_q) {
      N <- N + q[t] - prev_q
    }
    prev_q <- q[t]
  }
  
  # Run through the state sequence, and pull out the durations of all the
  # intervals for which the system is at or above a given state greater than 1.
  bursts <- data.frame(level=rep(NA, N), start=rep(offsets[1], N),
                       end=rep(offsets[1], N)) # Using offsets[1] for the type.
  burstcounter <- 0
  prev_q <- 0
  stack <- rep(NA, N) # Keeps track of which bursts are currently open.
  stackcounter <- 0
  for (t in 1:n) {
    if (q[t] > prev_q) {
      num_levels_opened <- q[t] - prev_q
      for (i in 1:num_levels_opened) {
        burstcounter <- burstcounter + 1
        bursts$level[burstcounter] <- prev_q + i
        bursts$start[burstcounter] <- offsets[t]
        stackcounter <- stackcounter + 1
        stack[stackcounter] <- burstcounter
      }
    } else if (q[t] < prev_q) {
      num_levels_closed <- prev_q - q[t]
      for (i in 1:num_levels_closed) {
        bursts$end[stack[stackcounter]] <- offsets[t]
        stackcounter <- stackcounter - 1
      }
    }
    prev_q <- q[t]
  }
  while (stackcounter >= 0) {
    bursts$end[stack[stackcounter]] <- offsets[n + 1]
    stackcounter <- stackcounter - 1
  }
  
  class(bursts) <- c('bursts', 'data.frame')
  bursts
}

plot.bursts <- function (x, ...)
{
  max_level <- max(x$level)
  
  plot(c(), type='n', xlab='Time', ylab='Level', bty='n',
       xlim=c(x$start[1], x$end[1]),
       ylim=c(1, max_level), yaxt='n', ...)
  axis(2, at=1:max_level)
  
  arrows(x$start, x$level, x$end, x$level,
         code=3, angle=90, length=0.05)
}
