# Monte Carlo Simulation
# Taken from "Eonometrics in R" by Grant V. Farnsworth
#
# The following block of code creates a randomly distributed data X with
# 25 members. It then creates a y vector that is conditionally distributed
# as
#       y = 2 + 3x + e
# It then does a regression of x on y and stores the slope coefficient.
# The generation of y and calculation of the slope coefficient are
# repeated 600 times. The mean and sample variance of the slope coefficient
# are then calculated. A comparison of the sample variance of the estimated
# coefficient with the analytic solution for the variance of the slope
# coefficient is then possible. 

init_func <- function() {
   set.seed(123) # comment this out if you don't want deterministic results
   x <<- rnorm(25, mean=2, sd=1)
}

comp_func <- function(...) {
   y <- rnorm(25, mean=(3*x+2), sd=1)
   beta <- lm(y~x)
   beta$coef[2]
}

if (TRUE) {
  cat('Parallel version\n')
  library(nws)
  s <- sleigh(workerCount=3)
  eachWorker(s, init_func)
  eo = list(chunkSize = 200)
  temp_results <- eachElem(s, comp_func, list(1:600), eo=eo) 
  A = unlist(temp_results)
  Abar <- mean(A)
  varA <- var(A)
  cat(Abar, varA, '\n')
}

if (TRUE) {
  cat('Sequential version\n')
  init_func()
  A <- vector(length=600)
  for (i in 1:600) {
    A[i] <- comp_func(i)
  }
  Abar <- mean(A)
  varA <- var(A)
  cat(Abar, varA, '\n')
}
