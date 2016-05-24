# test ODE solvers euler, midpoint and RK4 by
# solving dy/dt = -y, y(0) = 1, from t = 0 to t = 1
# we halve the step size until a given accuracy is achieved

rm(list = ls())
source("../scripts/euler.r")
source("../scripts/midpoint.r")
source("../scripts/RK4.r")

dydt <- function(t=NULL, y) -y  # y(t) = y(0)*exp(-t)

last <- function(x) x[length(x)]

# tol.vec is a vector of tolerances
# for each element of tol.vec a col of ftn.calls
# counts the function calls by each solver
tol.vec <- 10^(-(1:6))
ftn.calls <- matrix(nrow=3, ncol=length(tol.vec))
for (i in 1:length(tol.vec)) {
  tol <- tol.vec[i]
  # euler
  h <- 1
  y <- last(euler(dydt, 0, 1, h, n=1/h))
  while (abs(y - exp(-1)) > tol) {
    h <- h/2
    y <- last(euler(dydt, 0, 1, h, n=1/h))
  }
  ftn.calls[1,i] <- 1/h
  # midpoint
  h <- 1
  y <- last(midpoint(dydt, 0, 1, h, n=1/h))
  while (abs(y - exp(-1)) > tol) {
    h <- h/2
    y <- last(midpoint(dydt, 0, 1, h, n=1/h))
  }
  ftn.calls[2,i] <- 2/h
  # RK4
  h <- 1
  y <- last(RK4(dydt, 0, 1, h, n=1/h))
  while (abs(y - exp(-1)) > tol) {
    h <- h/2
    y <- last(RK4(dydt, 0, 1, h, n=1/h))
  }
  ftn.calls[3,i] <- 4/h
}

# print tolerance vs function calls for each solver
# on a log-log scale
plot(-log(tol.vec), log(ftn.calls[1,]), type='b', pch=1,
     xlab = '-log(tolerance)', ylab = 'log(function calls)')
lines(-log(tol.vec), log(ftn.calls[2,]), type='b', pch=2)
lines(-log(tol.vec), log(ftn.calls[3,]), type='b', pch=3)
legend('topleft', legend = c('Euler', 'Midpoint', 'RK4'),
       pch = c(1, 2, 3), inset = .1)
