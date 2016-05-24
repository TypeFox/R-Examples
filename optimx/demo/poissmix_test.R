options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)


###############################################################
cat("optimx test poissmix.loglik ...\n")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

poissmix.loglik <- function(p,y) {
i <- 0:(length(y)-1)
loglik <- y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) + 
		(1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
return ( -sum(loglik) )
}

####################

# Real data from Hasselblad (JASA 1969)
poissmix.dat <- data.frame(death=0:9, freq=c(162,267,271,185,111,61,27,8,3,1))

lo <- c(0.001,0,0)
hi <- c(0.999, Inf, Inf)

y <- poissmix.dat$freq
p <- runif(3,c(0.3,1,1),c(0.7,5,8))
system.time(ans.optx <- optimx(par=p, fn=poissmix.loglik, y=y, lower=lo, upper=hi,
    control=list(maxit=2500,save.failures=TRUE,all.methods=TRUE)))[1]

print(ans.optx)


