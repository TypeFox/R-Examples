# Here is yet another silly example that computes pi.
# It demonstrates the use of the "timeout" function to limit
# the length of time to compute the estimate.  It is also
# yet another demonstration of using vector operations within
# a foreach loop.

library(itertools)
library(foreach)

# Initialize variables
n <- 1000000  # length of vectors
t <- 60       # seconds to compute

# Create iterators of random numbers uniformly distributed over [-0.5, 0.5]
xit <- irunif(n=n, min=-0.5, max=0.5)
yit <- irunif(n=n, min=-0.5, max=0.5)

# Create a timeout iterator that just happens to return zeros
timer <- timeout(irepeat(0), time=t)

# Define a ".final" function that calculates pi from estimates of pi/4
calc.pi <- function(x) {
  cat(sprintf('computed %d estimates of pi/4\n', length(x)))
  4 * mean(x)
}

# foreach iterates over "timer" even though it is not named.
# It's not named because it is needed to break out of the
# loop, not for its value.
pi <- foreach(x=xit, y=yit, timer, .combine='c', .final=calc.pi) %do% {
  sum(sqrt(x*x + y*y) < 0.5) / n
}

cat(sprintf('Approximate value of pi after about %d seconds: %f\n', t, pi))
