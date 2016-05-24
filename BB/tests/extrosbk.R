if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1236)
old.seed <- setRNG(test.rng)
#iseed <- 1236  # this seed was used for tests conducted on March 25, 2008.  
#set.seed(iseed)

extrosbk <- function(x) {
n <- length(x)
f <- rep(NA, n)
j <- 2 * (1:(n/2))
jm1 <- j - 1
f[jm1] <- 10 * (x[j] - x[jm1]^2)
f[j] <-  1 - x[jm1]
f
}

p0 <- runif(50)
ans1 <- dfsane(par=p0, fn=extrosbk, method=1)
ans2 <- dfsane(par=p0, fn=extrosbk, method=2)
ans3 <- sane(par=p0, fn=extrosbk, method=1)
ans4 <- sane(par=p0, fn=extrosbk, method=2)
#ans <- nlsolve(p0, fn=extrosbk, method="L-BFGS-B")

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) #, ans$val) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) #, ans$counts[1]) 

nsim <- 10
dfsane1.extrosbk <- dfsane2.extrosbk <- sane1.extrosbk <- sane2.extrosbk <- matrix(NA, nsim, 5)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- runif(50)
t1 <- system.time(ans <- sane(par=p0, fn=extrosbk, method=1, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans))sane1.extrosbk[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t1)
t2 <- system.time(ans <- sane(par=p0, fn=extrosbk, method=2, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans))sane2.extrosbk[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t2)
t3 <- system.time(ans <- dfsane(par=p0, fn=extrosbk, method=1, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans))dfsane1.extrosbk[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t3)
t4 <- system.time(ans <- dfsane(par=p0, fn=extrosbk, method=2, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans)) dfsane2.extrosbk[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t4)
}

z <- apply(sane1.extrosbk, 2, summary)
print(z)
print(z[,1], digits=18)
z <- apply(sane2.extrosbk, 2, summary)
print(z)
print(z[,1], digits=18)
z <- apply(dfsane1.extrosbk, 2, summary)
print(z)
print(z[,1], digits=18)
z <- apply(dfsane2.extrosbk, 2, summary)
print(z)
print(z[,1], digits=18)

