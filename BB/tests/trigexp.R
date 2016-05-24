if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1236)
old.seed <- setRNG(test.rng)
#iseed <- 1236  # this seed was used for tests conducted on March 25, 2008.  
#set.seed(iseed)

trigexp <- function(x) {
n <- length(x)
F <- rep(NA, n)
F[1] <- 3*x[1]^2 + 2*x[2] - 5 + sin(x[1] - x[2]) * sin(x[1] + x[2])
tn1 <- 2:(n-1)
F[tn1] <- -x[tn1-1] * exp(x[tn1-1] - x[tn1]) + x[tn1] * ( 4 + 3*x[tn1]^2) +
        2 * x[tn1 + 1] + sin(x[tn1] - x[tn1 + 1]) * sin(x[tn1] + x[tn1 + 1]) - 8 
F[n] <- -x[n-1] * exp(x[n-1] - x[n]) + 4*x[n] - 3
F
}

p0 <- rnorm(50, sd=2)
ans1 <- dfsane(par=p0, fn=trigexp, method=1)
ans2 <- dfsane(par=p0, fn=trigexp, method=2)
ans3 <- sane(par=p0, fn=trigexp, method=2)
ans4 <- sane(par=p0, fn=trigexp, method=3)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 


nsim <- 50
dfsane1.trigexp <- dfsane2.trigexp <- sane1.trigexp <- sane2.trigexp <- matrix(NA, nsim, 5)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- rnorm(50)
t1 <- system.time(ans <- sane(par=p0, fn=trigexp, method=1, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans))sane1.trigexp[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t1)
t2 <- system.time(ans <- sane(par=p0, fn=trigexp, method=2, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans))sane2.trigexp[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t2)
t3 <- system.time(ans <- dfsane(par=p0, fn=trigexp, method=1, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans))dfsane1.trigexp[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t3)
t4 <- system.time(ans <- dfsane(par=p0, fn=trigexp, method=2, control=list(BFGS=TRUE, trace=F)))[1]
if (!is.null(ans)) dfsane2.trigexp[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t4)
}

z <- apply(sane1.trigexp, 2, summary)
print(z)
print(z[,1], digits=18)
z <- apply(sane2.trigexp, 2, summary)
print(z)
print(z[,1], digits=18)
z <- apply(dfsane1.trigexp, 2, summary)
print(z)
print(z[,1], digits=18)
z <- apply(dfsane2.trigexp, 2, summary)
print(z)
print(z[,1], digits=18)

