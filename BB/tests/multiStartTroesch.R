if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1234)
old.seed <- setRNG(test.rng)

troesch <- function(x) {
  n <- length(x)
  tnm1 <- 2:(n-1)
  f <- rep(NA, n)
    h <- 1 / (n+1)
    h2 <- 10 * h^2
    f[1] <- 2 * x[1] + h2 * sinh(10 * x[1]) - x[2] 
    f[tnm1] <- 2 * x[tnm1] + h2 * sinh(10 * x[tnm1]) - x[tnm1-1] - x[tnm1+1]    

    f[n] <- 2 * x[n] + h2 * sinh(10* x[n]) - x[n-1] - 1
  f
  }
  

p0 <- matrix(runif(50), 5, 10)  # 5 starting values, each of length 10
#ans <- BBsolve(par=p0, fn=troesch)
ans <- multiStart(par=p0, fn=troesch)
ans$par
ans$info

# ans$convergence
 
z <- sum(ans$par)
good   <-    2.191639378551349
#on Windows 
#on Linux64 
#on Linux32  2.191639378551349
print(z, digits=16)
if(any(abs(good - z) > 5e-9)) stop("BB test BBsolve Troesch FAILED")
