options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)


##########
cat("BB test dfsane ...\n")

expo1 <- function(x) {
    	n <- length(x)
    	f <- rep(NA, n)
    	f[1] <- exp(x[1] - 1) - 1
    	f[2:n] <- (2:n) * (exp(x[2:n] - 1) - x[2:n])
    	f
    	}

ans <- dfsane(par=runif(50), fn=expo1)
 
z <- sum(ans$par)
good   <-   49.9963301148
print(z, digits=16)
if(any(abs(good - z) > 1e-9)) stop("BB test dfsane FAILED")
