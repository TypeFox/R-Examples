
options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("BB test broydt.f ...\n")

broydt.f <- function(x) {
n <- length(x)
f <- rep(NA, n)
f[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
sum(f*f)
}

p0 <- rnorm(50, sd=1)
system.time(ans.spg <- spg(par=p0, fn=broydt.f, control=list(maxit=10000)))[1]
 
z <- sum(ans.spg$par)
# good   <-   13.6791297530 converged ?
good   <-   3.393768541175171
print(z, digits=16)
if(any(abs(good - z) > 5e-1)) stop("BB test broydt.f a FAILED")
 
