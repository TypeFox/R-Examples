options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)


##########
cat("BB test trig.f ...\n")

trig.f <- function(x){
n <- length(x)
i <- 1:n
f <- n - sum(cos(x)) + i*(1 - cos(x)) - sin(x) 
sum(f*f)
}

p0 <- rnorm(50,sd=5)
system.time(ans.spg <- spg(par=p0, fn=trig.f, control=list(maxit=2500)))[1]
 
z <- sum(ans.spg$par)
# -5.3788206334      Windows
# -5.379284954782386 Ubuntu-32 11.10
good   <-   -5.379
print(z, digits=16)
if(any(abs(good - z) > 1e-3)) stop("BB test trig.f a FAILED")
 
