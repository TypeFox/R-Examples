# J C Nash 2010-2-11 for optimx
# test results in these files are indicated for
# not yet included
# Ravi Varadhan version of Brown Almost Linear Function

options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx test brown-x.f ...\n")

brown.f <- function(x) {
p <- x
n <- length(p)
odd <- seq(1,n,by=2)
even <- seq(2,n,by=2)
sum((p[odd]^2)^(p[even]^2 + 1) + (p[even]^2)^(p[odd]^2 + 1))
}

brown.g <- function(x) { # gradient
p <- x
n <- length(p)
gg<-rep(0,n) # create the gradient vector
odd <- seq(1,n,by=2)
even <- seq(2,n,by=2)

sum((p[odd]^2)^(p[even]^2 + 1) + (p[even]^2)^(p[odd]^2 + 1))

  stop("NOT IMPLEMENTED YET")
}

npar<-50 # Down from 500
p0 <- rnorm(npar,sd=2)
system.time(ans.brownRV <- optimx(par=p0, fn=brown.f, 
   control=list(trace=2, all.methods=TRUE, save.failures=TRUE, maxit=2500)))[1]
print(ans.brownRV)

cat("================= end brownRV_test ==================\n")


