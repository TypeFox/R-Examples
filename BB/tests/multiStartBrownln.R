if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1234)
old.seed <- setRNG(test.rng)


# Brown's almost linear system(A.P. Morgan, ACM 1983)
brownlin <- function(x) {
# two distinct solutions if n is even
# three distinct solutions if n is odd
  	n <- length(x)
  	f <- rep(NA, n)
	nm1 <- 1:(n-1)
	f[nm1] <- x[nm1] + sum(x) - (n+1)
	f[n] <- prod(x) - 1 
	f
}

p0 <- matrix(rnorm(200), 20, 10)  # 20 starting values, each of length 10
#ans1 <- dfsane(par=p0[1,], fn=brownlin)
#ans <- BBsolve(par=p0, fn=brownlin)
ans <- multiStart(par=p0, fn=brownlin)

#ans$par

#pc <- princomp(ans$par) 
#plot(pc$scores[,1])  # plot shows two/three distinct solutions depending on n is even/odd


# ans$convergence
 
z <- sum(ans$par)
good   <-    200.2262649045772
#on Windows 
#on Linux64  200.205695196736
#on Linux32  200.2262649045772
print(z, digits=16)
if(any(abs(good - z) > 5e-1)) stop("BB test BBsolve Brownln FAILED")
