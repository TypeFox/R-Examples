options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

#########################################################################################
cat("BB test sc2 ...\n")

sc2.f <- function(x){
n <- length(x)
vec <- 1:n
sum(vec * (exp(x) - x)) / 10
}

sc2.g <- function(x){
n <- length(x)
vec <- 1:n
vec * (exp(x) - 1) / 10
}

neg.sc2.f <- function(x){
n <- length(x)
vec <- 1:n
-sum(vec * (exp(x) - x)) / 10
}

neg.sc2.g <- function(x){
n <- length(x)
vec <- 1:n
-vec * (exp(x) - 1) / 10
}

p0 <- runif(50, min=-1, max=1)
system.time(ans.spg <- spg(par=p0, fn=sc2.f, control=list(maxit=2500)))[1]

z <- sum(ans.spg$par)
# -0.0002158022390025393 on Windows 
#  0.001523050487259005 on Windows i386-pc-mingw32 (32-bit) R 2.13.0 (2011-04-13)

good   <-    0.0
print(z, digits=16)
if(any(abs(good - z) > 0.002)) stop("BB test sc2 a FAILED")

system.time(neg.ans.spg <- spg(par=p0, fn=neg.sc2.f, 
              control=list(maxit=2500, maximize=TRUE)))[1]

z <- sum(neg.ans.spg$par)
good   <-    0.0
print(z, digits=16)
if(any(abs(good - z) > 0.002)) stop("BB test neg sc2 a FAILED")

system.time(ans.spg <- spg(par=p0, fn=sc2.f, gr=sc2.g,
   control=list(maxit=2500)))[1]

z <- sum(ans.spg$par)
# 2.565413040899874e-06 Linux64 (mfacl2)
# 6.677493403589264e-05 Linux64
# 0.0002100097368926836 i386-pc-solaris2.10 (32-bit) CRAN R 2.13.0 Patched (2011-04-15 r55454)
# 0.0006271870150276102 i386-apple-darwin9.8.0 (32-bit) CRAN R 2.13.0 Under development (unstable) (2011-03-07 r54691) 
# 0.0005317982091112333 x86_64-unknown-linux-gnu (64-bit) CRAN fedora 2.14.0 Under development (unstable) (2011-04-14 r55450) 

good <- 0.0
print(z, digits=16)
# test tol relaxed from 1e-4 to 1e-3 when ftol arg added to spg  2011.2-1
if(any(abs(good - z) >  0.001)) stop("BB test sc2 b FAILED")

system.time(neg.ans.spg <- spg(par=p0, fn=neg.sc2.f, gr=neg.sc2.g,
   control=list(maxit=2500, maximize=TRUE)))[1]

z <- sum(neg.ans.spg$par)
good <- 0.0
print(z, digits=16)
# test tol relaxed from 1e-4 to 1e-3 when ftol arg added to spg  2011.2-1
if(any(abs(good - z) >  0.001)) stop("BB test neg.sc2 b FAILED")

