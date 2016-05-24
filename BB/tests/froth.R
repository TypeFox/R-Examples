options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)


#########################################
cat("BB test froth ...\n")

froth <- function(p){
# Freudenstein and Roth function (Broyden, Mathematics of Computation 1965, p. 577-593)
f <- rep(NA,length(p))
f[1] <- -13 + p[1] + (p[2]*(5 - p[2]) - 2) * p[2]
f[2] <- -29 + p[1] + (p[2]*(1 + p[2]) - 14) * p[2]
sum (f * f)
}

p0 <- rpois(2,10)
system.time(ans.spg <- spg(par=p0, fn=froth, control=list(M=20, maxit=2500)))[1]
#ans.spg
system.time(ans.opt <- optim(par=p0, fn=froth, method="L-BFGS-B"))[1]

z <- sum(ans.spg$par)
good   <-   10.51597043896899
#on Windows 10.51597043896899
#on Linux64 10.51597043896899
#on Linux64 10.51597048420193 Ubuntu 10.04
#on Linux32 10.51597043896899
print(z, digits=16)
if(any(abs(good - z) > 1e-7)) stop("BB test froth a FAILED")
 
z <- sum(ans.opt$par)
good    <-  9.00000701456296
#on Windows 9.00000701456296
#on Linux64 9.00000701456296
#on Linux32 9.00000701456296
print(z, digits=16)
if(any(abs(good - z) > 1e-12)) stop("BB test froth  b FAILED")
