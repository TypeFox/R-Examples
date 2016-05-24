# J C Nash 2010-2-11 for optimx
# test results in these files are indicated for
# not yet included

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

npar<-50 # Down from 500
p0 <- rnorm(npar,sd=2)
system.time(ans.optx <- optimx(par=p0, fn=brown.f, control=list(all.methods=TRUE, save.failures=TRUE, maxit=2500)))[1]


print(ans.optx)

#allpar<-ans.optx$par # ans.optx is a dataframe!
#allmeth<-ans.optx$method
#nanswer<-length(allpar)
#
#for (i in 1:nanswer) {
#	curmeth<-allmeth[[i]]
#	z <- sum(ans.optx$par[[i]])
#	cat(curmeth,": ")
#	print(z, digits=16)
#}
 

