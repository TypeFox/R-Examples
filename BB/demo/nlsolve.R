if(!require("BB"))    stop("this requires package BB.")
if(!require("setRNG"))stop("this requires setRNG.")

#nsim <- 100
# nsim reduced from 100 to 20 because of testing time constraints on CRAN
nsim <- 20

# This was used for tests conducted on March 25, 2008, using setRNG(test.rng).  
#   iseed <- 1236  
# Replaced with setRNG to ensure rng and normal generators are set too.
# Changed from kind="Wichmann-Hill", normal.kind="Box-Muller", 
#  (back) to "Mersenne-Twister", normal.kind="Inversion", Jan 15, 2009
test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1236)
old.seed <- setRNG(test.rng)

#
# Some examples illustrating the use of SANE & DFSANE
#
expo1 <- function(x) {
#  From La Cruz and Raydan, Optim Methods and Software 2003, 18 (583-599)
n <- length(x)
f <- rep(NA, n)
f[1] <- exp(x[1] - 1) - 1
f[2:n] <- (2:n) * (exp(x[2:n] - 1) - x[2:n])
f
}
p0 <- rnorm(100)
ans1 <- dfsane(par=p0, fn=expo1, method=1)
ans2 <- dfsane(par=p0, fn=expo1, method=2)
ans3 <- dfsane(par=p0, fn=expo1, method=3)
ans4 <- sane(par=p0, fn=expo1)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 

setRNG(test.rng)
dfsane1.expo1 <- dfsane2.expo1 <- dfsane3.expo1 <- sane.expo1 <- matrix(NA, nsim, 4)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- rnorm(100)
ans <- sane(par=p0, fn=expo1, control=list(trace=F))
if (!is.null(ans)) sane.expo1[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=expo1, method=1, control=list(trace=F))
if (!is.null(ans)) dfsane1.expo1[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=expo1, method=2, control=list(trace=F))
if (!is.null(ans)) dfsane2.expo1[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=expo1, method=3, control=list(trace=F))
if (!is.null(ans)) dfsane3.expo1[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
}

apply(sane.expo1, 2, summary)
apply(dfsane1.expo1, 2, summary)
apply(dfsane2.expo1, 2, summary)
apply(dfsane3.expo1, 2, summary)

#expo1.results <- list(dfsane1=dfsane1.expo1, dfsane2=dfsane2.expo1, dfsane3=dfsane3.expo1, sane=sane.expo1) 
#dput(expo1.results, file="e:/bb/package/expo1.results")
# expo1.results <- dget(file="e:/bb/package/expo1.results")

############
expo3 <- function(p) {
#  From La Cruz and Raydan, Optim Methods and Software 2003, 18 (583-599)
n <- length(p)
f <- rep(NA, n)
onm1 <- 1:(n-1) 
f[onm1] <- onm1/10 * (1 - p[onm1]^2 - exp(-p[onm1]^2))
f[n] <- n/10 * (1 - exp(-p[n]^2))
f
}

n <- 100
p0 <- (1:n)/(4*n^2)
p0 <- rnorm(n, sd=4)
ans1 <- dfsane(par=p0, fn=expo3, method=1)
ans2 <- dfsane(par=p0, fn=expo3, method=2)
ans3 <- dfsane(par=p0, fn=expo3, method=3)
ans4 <- sane(par=p0, fn=expo3)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 

setRNG(test.rng)

dfsane1.expo3 <- dfsane2.expo3 <- dfsane3.expo3 <- sane.expo3 <- matrix(NA, nsim, 4)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- rnorm(100)
ans <- sane(par=p0, fn=expo3, control=list(trace=F))
if (!is.null(ans)) sane.expo3[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=expo3, method=1, control=list(trace=F))
if (!is.null(ans)) dfsane1.expo3[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=expo3, method=2, control=list(trace=F))
if (!is.null(ans)) dfsane2.expo3[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=expo3, method=3, control=list(trace=F))
if (!is.null(ans)) dfsane3.expo3[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
}

apply(sane.expo3, 2, summary)
apply(dfsane1.expo3, 2, summary)
apply(dfsane2.expo3, 2, summary)
apply(dfsane3.expo3, 2, summary)

#expo3.results <- list(dfsane1=dfsane1.expo3, dfsane2=dfsane2.expo3, dfsane3=dfsane3.expo3, sane=sane.expo3) 
#dput(expo3.results, file="e:/bb/package/expo3.results")

####################################################
froth <- function(p){
# Freudenstein and Roth function (Broyden, Mathematics of Computation 1965, p. 577-593)
f <- rep(NA,length(p))
f[1] <- -13 + p[1] + (p[2]*(5 - p[2]) - 2) * p[2]
f[2] <- -29 + p[1] + (p[2]*(1 + p[2]) - 14) * p[2]
f
}

p0 <- c(3,2)  # this gives the zero of the system
ans1 <- dfsane(par=p0, fn=froth, method=1)
ans2 <- dfsane(par=p0, fn=froth, method=2)
ans3 <- dfsane(par=p0, fn=froth, method=3)
ans4 <- sane(par=p0, fn=froth)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 


p0 <- c(1,1)  # this gives the local minimum that is not the zero of the system
ans1 <- dfsane(par=p0, fn=froth, method=1)
ans2 <- dfsane(par=p0, fn=froth, method=2)
ans3 <- dfsane(par=p0, fn=froth, method=3)
ans4 <- sane(par=p0, fn=froth)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 


p0 <- rpois(2,10)  # trying random starts
ans1 <- dfsane(par=p0, fn=froth, method=1)
ans2 <- dfsane(par=p0, fn=froth, method=2)
ans3 <- dfsane(par=p0, fn=froth, method=3)
ans4 <- sane(par=p0, fn=froth)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 

###########################################
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

p0 <- rnorm(100, sd=3)
ans1 <- dfsane(par=p0, fn=trigexp, method=1)
ans2 <- dfsane(par=p0, fn=trigexp, method=2)
ans3 <- dfsane(par=p0, fn=trigexp, method=3)
ans4 <- sane(par=p0, fn=trigexp)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 


setRNG(test.rng)

dfsane1.trigexp <- dfsane2.trigexp <- dfsane3.trigexp <- sane.trigexp <- matrix(NA, nsim, 4)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- rnorm(100)
ans <- sane(par=p0, fn=trigexp)
if (!is.null(ans)) sane.trigexp[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=trigexp, method=1)
if (!is.null(ans)) dfsane1.trigexp[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=trigexp, method=2)
if (!is.null(ans)) dfsane2.trigexp[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=trigexp, method=3)
if (!is.null(ans)) dfsane3.trigexp[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
}

apply(sane.trigexp, 2, summary)
apply(dfsane1.trigexp, 2, summary)
apply(dfsane2.trigexp, 2, summary)
apply(dfsane3.trigexp, 2, summary)

#trigexp.results <- list(dfsane1=dfsane1.trigexp, dfsane2=dfsane2.trigexp, dfsane3=dfsane3.trigexp, sane=sane.trigexp) 
#dput(trigexp.results, file="e:/bb/package/trigexp.results")


###########################################
valley <- function(x) {
c1 <- 1.003344481605351
c2 <- -3.344481605351171e-04
n <- length(x)
f <- rep(NA, n)
j <- 3 * (1:(n/3))
jm2 <- j - 2
jm1 <- j - 1
f[jm2] <- (c2 * x[jm2]^3 + c1 * x[jm2]) * exp(-(x[jm2]^2)/100) - 1
f[jm1] <- 10 * (sin(x[jm2]) - x[jm1])
f[j] <- 10 * (cos(x[jm2]) - x[j])
f
}
 
p0 <- rnorm(102, sd=3)  # number of unknowns must be a multiple of 3
ans1 <- dfsane(par=p0, fn=valley, method=1)
ans2 <- dfsane(par=p0, fn=valley, method=2)
ans3 <- dfsane(par=p0, fn=valley, method=3)
ans4 <- sane(par=p0, fn=valley)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 


setRNG(test.rng)

dfsane1.valley <- dfsane2.valley <- dfsane3.valley <- sane.valley <- matrix(NA, nsim, 4)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- rnorm(102)
ans <- sane(par=p0, fn=valley, control=list(trace=F))
if (!is.null(ans)) sane.valley[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=valley, method=1, control=list(trace=F))
if (!is.null(ans)) dfsane1.valley[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=valley, method=2, control=list(trace=F))
if (!is.null(ans)) dfsane2.valley[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=valley, method=3, control=list(trace=F))
if (!is.null(ans)) dfsane3.valley[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
}

apply(sane.valley, 2, summary)
apply(dfsane1.valley, 2, summary)
apply(dfsane2.valley, 2, summary)
apply(dfsane3.valley, 2, summary)

#valley.results <- list(dfsane1=dfsane1.valley, dfsane2=dfsane2.valley, dfsane3=dfsane3.valley, sane=sane.valley) 
#dput(valley.results, file="e:/bb/package/valley.results")


###########################################
broydt <- function(x) {
n <- length(x)
f <- rep(NA, n)
f[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
f
}

p0 <- rnorm(500, sd=5)
ans1 <- dfsane(par=p0, fn=broydt, method=1)
ans2 <- dfsane(par=p0, fn=broydt, method=2)
ans3 <- dfsane(par=p0, fn=broydt, method=3)
ans4 <- sane(par=p0, fn=broydt)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 

setRNG(test.rng)

dfsane1.broydt <- dfsane2.broydt <- dfsane3.broydt <- sane.broydt <- matrix(NA, nsim, 4)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- rnorm(100)
ans <- sane(par=p0, fn=broydt, control=list(trace=F))
if (!is.null(ans)) sane.broydt[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=broydt, method=1, control=list(trace=F))
if (!is.null(ans)) dfsane1.broydt[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=broydt, method=2, control=list(trace=F))
if (!is.null(ans)) dfsane2.broydt[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=broydt, method=3, control=list(trace=F))
if (!is.null(ans)) dfsane3.broydt[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
}

apply(sane.broydt, 2, summary)
apply(dfsane1.broydt, 2, summary)
apply(dfsane2.broydt, 2, summary)
apply(dfsane3.broydt, 2, summary)

#broydt.results <- list(dfsane1=dfsane1.broydt, dfsane2=dfsane2.broydt, dfsane3=dfsane3.broydt, sane=sane.broydt) 
#dput(broydt.results, file="e:/bb/package/broydt.results")

######################################
brent <- function(x) {
  n <- length(x)
  tnm1 <- 2:(n-1)
  F <- rep(NA, n)

	F[1] <- 3 * x[1] * (x[2] - 2*x[1]) + (x[2]^2)/4 
	F[tnm1] <-  3 * x[tnm1] * (x[tnm1+1] - 2 * x[tnm1] + x[tnm1-1]) + ((x[tnm1+1] - x[tnm1-1])^2) / 4	
	F[n] <- 3 * x[n] * (20 - 2 * x[n] + x[n-1]) +  ((20 - x[n-1])^2) / 4
  F
  }

p0 <- rnorm(100)
ans1 <- dfsane(par=p0, fn=brent, method=1)
ans2 <- dfsane(par=p0, fn=brent, method=2)
ans3 <- dfsane(par=p0, fn=brent, method=3)
ans4 <- sane(par=p0, fn=brent)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 


setRNG(test.rng)

dfsane1.brent <- dfsane2.brent <- dfsane3.brent <- sane.brent <- matrix(NA, nsim, 4)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
p0 <- rnorm(100)
ans <- sane(par=p0, fn=brent, control=list(trace=F))
if (!is.null(ans))sane.brent[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=brent, method=1, control=list(trace=F))
if (!is.null(ans))dfsane1.brent[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=brent, method=2, control=list(trace=F))
if (!is.null(ans)) dfsane2.brent[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=brent, method=3, control=list(trace=F))
if (!is.null(ans)) dfsane3.brent[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
}

apply(sane.brent, 2, summary)
apply(dfsane1.brent, 2, summary)
apply(dfsane2.brent, 2, summary)
apply(dfsane3.brent, 2, summary)

#brent.results <- list(dfsane1=dfsane1.brent, dfsane2=dfsane2.brent, dfsane3=dfsane3.brent, sane=sane.brent) 
#dput(brent.results, file="e:/bb/package/brent.results")

######################################
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
  
p0 <- rnorm(100)
ans1 <- dfsane(par=p0, fn=troesch, method=1)
ans2 <- dfsane(par=p0, fn=troesch, method=2)
ans3 <- dfsane(par=p0, fn=troesch, method=3)
ans4 <- sane(par=p0, fn=troesch)

c(ans1$resid, ans2$resid,ans3$resid, ans4$resid) 
c(ans1$feval, ans2$feval,ans3$feval,ans4$feval) 


setRNG(test.rng)
dfsane1.troesch <- dfsane2.troesch <- dfsane3.troesch <- sane.troesch <- matrix(NA, nsim, 4)
for (i in 1:nsim) {
cat("Simulation" , i, "\n")
#p0 <- rnorm(100)  # this doesn't work for "sane"; but works well for "dfsane"  
p0 <- runif(100)  # this works for all schemes
ans <- sane(par=p0, fn=troesch, control=list(trace=F))
if (!is.null(ans)) sane.troesch[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=troesch, method=1, control=list(trace=F))
if (!is.null(ans)) dfsane1.troesch[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=troesch, method=2, control=list(trace=F))
if (!is.null(ans)) dfsane2.troesch[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
ans <- dfsane(par=p0, fn=troesch, method=3, control=list(trace=F))
if (!is.null(ans)) dfsane3.troesch[i, ] <- c(ans$resid, ans$feval, ans$iter, ans$conv)
}

apply(sane.troesch, 2, summary)
apply(dfsane1.troesch, 2, summary)
apply(dfsane2.troesch, 2, summary)
apply(dfsane3.troesch, 2, summary)

#troesch.results <- list(dfsane1=dfsane1.troesch, dfsane2=dfsane2.troesch, dfsane3=dfsane3.troesch, sane=sane.troesch) 
#dput(troesch.results, file="e:/bb/package/troesch.results")
# troesch.results <- dget(file="e:/bb/package/troesch.results")

######################################
