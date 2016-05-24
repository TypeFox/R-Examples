# commented out optim examples because of time limitations testing on CRAN

###################################################
if(!require("BB"))    stop("this requires package BB.")
if(!require(numDeriv))stop("this requires package numDeriv.")
if(!require("setRNG"))stop("this requires setRNG.")

# This was used for tests conducted on March 25, 2008, using set.seed(test.rng).  
#   iseed <- 1236  
# Replaced April 7, 2008, with setRNG to ensure rng and normal generators are set too.
# Changed from kind="Wichmann-Hill", normal.kind="Box-Muller", 
#  (back) to "Mersenne-Twister", normal.kind="Inversion", Jan 15, 2009

test.rng <- list(kind="Mersenne-Twister", normal.kind="Inversion", seed=1236)
old.seed <- setRNG(test.rng)


fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}

rosbkext.f <- function(x){
p <- x
n <- length(p)
sum (100*(p[1:(n-1)]^2 - p[2:n])^2 + (p[1:(n-1)] - 1)^2)
}

sc2.f <- function(x){
nvec <- 1:length(x)
sum(nvec * (exp(x) - x)) / 10
}

sc2.g <- function(x){
nvec <- 1:length(x)
nvec * (exp(x) - 1) / 10
}

trig.f <- function(x){
n <- length(x)
i <- 1:n
f <- n - sum(cos(x)) + i*(1 - cos(x)) - sin(x) 
sum(f*f)
}

brown.f <- function(x) {
p <- x
n <- length(p)
odd <- seq(1,n,by=2)
even <- seq(2,n,by=2)
sum((p[odd]^2)^(p[even]^2 + 1) + (p[even]^2)^(p[odd]^2 + 1))
}

froth <- function(p){
# Freudenstein and Roth function (Broyden, Mathematics of Computation 1965, p. 577-593)
f <- rep(NA,length(p))
f[1] <- -13 + p[1] + (p[2]*(5 - p[2]) - 2) * p[2]
f[2] <- -29 + p[1] + (p[2]*(1 + p[2]) - 14) * p[2]
sum (f * f)
}

chen.f <- function(x) {
v <- log(x) + exp(x)
f <- (v - sqrt(v^2 + 5e-04))/2
sum (f * f)
}

valley.f <- function(x) {
c1 <- 1.003344481605351
c2 <- -3.344481605351171e-03
n <- length(x)
f <- rep(NA, n)
j <- 3 * (1:(n/3))
jm2 <- j - 2
jm1 <- j - 1
f[jm2] <- (c2 * x[jm2]^3 + c1 * x[jm2]) * exp(-(x[jm2]^2)/100) - 1
f[jm1] <- 10 * (sin(x[jm2]) - x[jm1])
f[j] <- 10 * (cos(x[jm2]) - x[j])
sum(f*f)
}
 
broydt.f <- function(x) {
n <- length(x)
f <- rep(NA, n)
f[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
sum(f*f)
}

#########################################################################################
p0 <- rnorm(2,sd=2)
system.time(ans.spg <- spg(par=p0, fn=fr))[1]
system.time(ans.spg <- spg(par=p0, fn=fr, method=1))[1]
system.time(ans.spg <- spg(par=p0, fn=fr, method=2))[1]
#system.time(ans.opt <- optim(par=p0, fn=fr, method="L-BFGS-B"))[1]
#########################################################################################
p0 <- rnorm(200,sd=2)
system.time(ans.spg <- spg(par=p0, fn=sc2.f))[1]
system.time(ans.spg <- spg(par=p0, fn=sc2.f, method=1))[1]
system.time(ans.spg <- spg(par=p0, fn=sc2.f, method=2))[1]
#system.time(ans.opt <- optim(par=p0, fn=sc2.f, method="L-BFGS-B"))[1]

# This demonstrates the value of providing "exact" gradient information
# Much faster computation, as well as better convergence
system.time(ans.spg <- spg(par=p0, fn=sc2.f, gr=sc2.g))[1]  

##########
p0 <- rnorm(200,sd=2)
system.time(ans.spg <- spg(par=p0, fn=brown.f))[1]
system.time(ans.spg <- spg(par=p0, fn=brown.f, meth=1))[1]
system.time(ans.spg <- spg(par=p0, fn=brown.f, method=2))[1]
#system.time(ans.opt <- optim(par=p0, fn=brown.f, method="L-BFGS-B"))[1]

##########
p0 <- rnorm(200,sd=2)
system.time(ans.spg <- spg(par=p0, fn=rosbkext.f))[1]
#system.time(ans.spg <- spg(par=p0, fn=rosbkext.f, meth=1))[1]
#system.time(ans.spg <- spg(par=p0, fn=rosbkext.f, method=2))[1]
#system.time(ans.opt <- optim(par=p0, fn=rosbkext.f, method="L-BFGS-B"))[1]
 
##########
p0 <- rnorm(200,sd=5)
system.time(ans.spg <- spg(par=p0, fn=trig.f))[1]
system.time(ans.spg <- spg(par=p0, fn=trig.f, meth=1))[1]
system.time(ans.spg <- spg(par=p0, fn=trig.f, method=2))[1]
#system.time(ans.opt <- optim(par=p0, fn=trig.f, method="L-BFGS-B"))[1]

##########
p0 <- rexp(500)
system.time(ans.spg <- spg(par=p0, fn=chen.f, lower=0))[1]
system.time(ans.spg <- spg(par=p0, fn=chen.f, lower=0, meth=1))[1]
system.time(ans.spg <- spg(par=p0, fn=chen.f, lower=0, method=2))[1]
#system.time(ans.opt <- optim(par=p0, fn=chen.f, lower=0, method="L-BFGS-B"))[1]
##########
p0 <- rnorm(99, sd=2)
system.time(ans.spg <- spg(par=p0, fn=valley.f))[1]
system.time(ans.spg <- spg(par=p0, fn=valley.f, meth=1))[1]
system.time(ans.spg <- spg(par=p0, fn=valley.f, method=2))[1]
#system.time(ans.opt <- optim(par=p0, fn=valley.f, method="L-BFGS-B"))[1]

##########
p0 <- rnorm(500, sd=2)
system.time(ans.spg <- spg(par=p0, fn=broydt.f))[1]
system.time(ans.spg <- spg(par=p0, fn=broydt.f, meth=1))[1]
system.time(ans.spg <- spg(par=p0, fn=broydt.f, method=2))[1]
#system.time(ans.opt <- optim(par=p0, fn=broydt.f, method="L-BFGS-B"))[1]
#########################################
p0 <- rpois(2,10)
ans.spg <- spg(par=p0, fn=froth)
ans.spg
ans.spg <- spg(par=p0, fn=froth, meth=1)
ans.spg
ans.spg <- spg(par=p0, fn=froth, method=2)
ans.spg
#optim(par=p0, fn=froth, method="L-BFGS-B")

###############################################################
poissmix.loglik <- function(p,y) {
i <- 0:(length(y)-1)
loglik <- y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) + 
		(1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
return (sum(loglik) )
}

###############################################################
# Real data from Hasselblad (JASA 1969)
poissmix.dat <- data.frame(death=0:9, freq=c(162,267,271,185,111,61,27,8,3,1))

lo <- c(0.001,0,0)
hi <- c(0.999, Inf, Inf)

y <- poissmix.dat$freq

p0 <- runif(3,c(0.2,1,1),c(0.8,5,8))  # randomly generated starting values
t.spg <- system.time(ans.spg <- spg(par=p0, fn=poissmix.loglik, y=y, 
   projectArgs=list(lower=lo, upper=hi), control=list(maximize=T)))[1]

t.spg <- system.time(ans.spg <- spg(par=p0, fn=poissmix.loglik, y=y, 
   projectArgs=list(lower=lo, upper=hi), control=list(maximize=T), meth=1))[1]

t.spg <- system.time(ans.spg <- spg(par=p0, fn=poissmix.loglik, y=y, 
   projectArgs=list(lower=lo, upper=hi), control=list(maximize=T), meth=2))[1]

#ans.opt <- optim(par=p0, fn=poissmix.loglik, y=y, method="L-BFGS-B", lower=lo, upper=hi, 
#	control=list(fnscale=-1))

grad(ans.spg$par, func=poissmix.loglik, y=y)
#grad(ans.opt$par, func=poissmix.loglik, y=y)

###############################################################
##
dvm <- function (theta, mu, kappa) 
{
    1/(2 * pi * besselI(x = kappa, nu = 0, expon.scaled = TRUE)) * 
        (exp(cos(theta - mu) - 1))^kappa
}

##
rmixedvm <- function (n, mu1, mu2, kappa1, kappa2, p) {
temp <- runif(n)
n1 <- sum(temp <= p)
y <- c(rvm(n1,mu1,kappa1),rvm(n-n1,mu2,kappa2))
return(y)
}

##
rvm <- function (n, mean, k) 
{
    vm <- c(1:n)
    a <- 1 + (1 + 4 * (k^2))^0.5
    b <- (a - (2 * a)^0.5)/(2 * k)
    r <- (1 + b^2)/(2 * b)
    obs <- 1
    while (obs <= n) {
        U1 <- runif(1, 0, 1)
        z <- cos(pi * U1)
        f <- (1 + r * z)/(r + z)
        c <- k * (r - f)
        U2 <- runif(1, 0, 1)
        if (c * (2 - c) - U2 > 0) {
            U3 <- runif(1, 0, 1)
            vm[obs] <- sign(U3 - 0.5) * acos(f) + mean
            vm[obs] <- vm[obs]%%(2 * pi)
            obs <- obs + 1
        }
        else {
            if (log(c/U2) + 1 - c >= 0) {
                U3 <- runif(1, 0, 1)
                vm[obs] <- sign(U3 - 0.5) * acos(f) + mean
                vm[obs] <- vm[obs]%%(2 * pi)
                obs <- obs + 1
            }
        }
    }
    vm
}

#
vmmix.loglik <- function(x,y){
p <- x
sum(log(p[5]*dvm(y,p[1],p[2])+(1-p[5])*dvm(y,p[3],p[4])))
}

y <- rmixedvm(n=500, mu1=pi/2, mu2=3*pi/2, kappa1=1.9, kappa2=2.2, p=0.67)
p <- c(pi/4,2,pi,1,0.5)

lo <- rep(0.001,5)
hi <- c(Inf, Inf, Inf, Inf, 0.999)

p0 <- c(runif(5,c(0,0.1,0,0.1,0.2),c(2*pi,5,2*pi,5,0.8)))

t.spg <- system.time(ans.spg <- spg(par=p0, fn=vmmix.loglik, y=y, 
   projectArgs=list(lower=lo, upper=hi), control=list(maximize=T, M=20)))[1]

t.spg <- system.time(ans.spg <- spg(par=p0, fn=vmmix.loglik, y=y, 
   projectArgs=list(lower=lo, upper=hi), method=1, 
   control=list(maximize=T, M=20)))[1]

t.spg <- system.time(ans.spg <- spg(par=p0, fn=vmmix.loglik, y=y, 
   projectArgs=list(lower=lo, upper=hi), method=2, 
   control=list(maximize=T, M=20)))[1]

#ans.opt <- optim(par=p0, fn=vmmix.loglik, y=y, method="L-BFGS-B", lower=lo, upper=hi, 
#  control=list(fnscale=-1))

grad(ans.spg$par, func=vmmix.loglik, y=y)
#grad(ans.opt$par, func=vmmix.loglik, y=y)

