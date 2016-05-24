if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1234)
old.seed <- setRNG(test.rng)


broydt <- function(x) {
n <- length(x)
f <- rep(NA, n)
#h <- 0.5
h <- 2
f[1] <- ((3 - h*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - h*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - h*x[n]) * x[n]) - x[n-1] + 1
f
}


p0 <- rnorm(50)
ans.opt <- BBsolve(par=p0, fn=broydt)  # note that the default doesn't work generally.

ans.opt$convergence
 
z <- sum(ans.opt$par)
#good   <-  -352.9190956497645 
good   <-  -34.72104 
#on Windows 
#on Linux64 
#on Linux32  -352.9190956497645
print(z, digits=16)
if(any(abs(good - z) > 5e-1)) stop("BB test BBsolve FAILED")
