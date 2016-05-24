options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)


##########
cat("BB test valley.f ...\n")

valley.f <- function(x) {
c1 <- 1.003344481605351
c2 <- -3.344481605351171e-03
n <- length(x)
f <- rep(NA, n)
j <- 3 * (1:(n/3))
jm2 <- j - 2
jm1 <- j - 1
f[jm2] <- (c2 * x[jm2]^3 + c1 * x[jm2]) * exp(-(x[jm2]^2)/100) - 1
f[jm1] <- 10 * (sin(x[jm2]) - x[jm1])
f[j] <- 10 * (cos(x[jm2]) - x[j])
sum(f*f)
}

p0 <- rnorm(99, sd=1)
system.time(ans.spg <- spg(par=p0, fn=valley.f))[1]
 
z <- sum(ans.spg$par)
good   <-   78.8339837886
print(z, digits=16)
if(any(abs(good - z) > 1e-2)) stop("BB test valley.f b FAILED")
 
