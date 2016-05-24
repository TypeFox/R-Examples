if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1236)
old.seed <- setRNG(test.rng)
#iseed <- 1236  # this seed was used for tests conducted on March 25, 2008.  
#set.seed(iseed)

troesch <- function(x) {
  n <- length(x)
  tnm1 <- 2:(n-1)
  F <- rep(NA, n)
    h <- 1 / (n+1)
    h2 <- 10 * h^2
    F[1] <- 2 * x[1] + h2 * sinh(10 * x[1]) - x[2] 
    F[tnm1] <- 2 * x[tnm1] + h2 * sinh(10 * x[tnm1]) - x[tnm1-1] - x[tnm1+1]    

    F[n] <- 2 * x[n] + h2 * sinh(10* x[n]) - x[n-1] - 1
  F
  }
  
p0 <- sort(runif(100))
ans1 <- dfsane(par=p0, fn=troesch, method=1)
ans2 <- dfsane(par=p0, fn=troesch, method=2)
ans3 <- sane(par=p0, fn=troesch, method=1)
ans4 <- sane(par=p0, fn=troesch, method=2)
#ans <- nlsolve(par=p0, fn=troesch)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 

# switched BFGS=TRUE to BFGS=FALSE below for speed, and 500 to 100

nsim <- 50
dfsane1.troesch <- dfsane2.troesch <- sane1.troesch <- sane2.troesch <- matrix(NA, nsim, 5)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- sort(runif(50))
t1 <- system.time(ans <- sane(par=p0, fn=troesch, method=1,
               control=list(BFGS=FALSE, trace=F)))[1]
if (!is.null(ans))sane1.troesch[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t1)
t2 <- system.time(ans <- sane(par=p0, fn=troesch, method=2,
                control=list(BFGS=FALSE, trace=F)))[1]
if (!is.null(ans))sane2.troesch[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t2)
t3 <- system.time(ans <- dfsane(par=p0, fn=troesch, method=1,
                control=list(BFGS=FALSE, trace=F)))[1]
if (!is.null(ans))dfsane1.troesch[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t3)
t4 <- system.time(ans <- dfsane(par=p0, fn=troesch, method=2,
                control=list(BFGS=FALSE, trace=F)))[1]
if (!is.null(ans)) dfsane2.troesch[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv, t4)
}

z <- apply(sane1.troesch, 2, summary)
print(z)
print(z[,1], digits=18)

z <- apply(sane2.troesch, 2, summary)
print(z)
print(z[,1], digits=18)

z <- apply(dfsane1.troesch, 2, summary)
print(z)
print(z[,1], digits=18)

z <- apply(dfsane2.troesch, 2, summary)
print(z)
print(z[,1], digits=18)


