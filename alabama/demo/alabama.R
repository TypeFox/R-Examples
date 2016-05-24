
#######################################################################################################
fn <- function(x) sum(x^2)

gr <- function(x) 2 * x

heq <- function(x) x[1]^2 - x[2] - x[3]^2 - 1

heq.jac <- function(x) {
j <- matrix(NA, 1, length(x))
j[1, 1] <- 2 * x[1]
j[1, 2] <- -1
j[1, 3] <- -2 * x[3]
j
}

set.seed(21) 
system.time(ans <- constrOptim.nl(par=rnorm(3), fn=fn, gr=gr, heq=heq, heq.jac=heq.jac))[1] 
ans

system.time(ans2 <- auglag(par=rnorm(3), fn=fn, gr=gr, heq=heq, heq.jac=heq.jac))[1] 
ans2
#######################################################################################################
fn <- function(x) (x[1] - 3/2)^2 + (x[2] - 1/8)^4

gr <- function(x) c(2 * (x[1] - 3/2) , 4 * (x[2] - 1/8)^3)

hin <- function(x) {
h <- rep(NA, 4)
h[1] <- 1 - x[1] - x[2]
h[2] <- 1 - x[1] + x[2]
h[3] <- 1 + x[1] - x[2]
h[4] <- 1 + x[1] + x[2]
h
}


hin.jac <- function(x) {
j <- matrix(NA, 4, length(x))
j[1, ] <- c(-1, -1)
j[2, ] <- c(-1, 1)
j[3, ] <- c(1, -1)
j[4, ] <- c(1, 1)
j
}

p0 <- runif(2) # note: initial value must be feasible
#p0 <- c(0.2, 0.6)
ans1 <- constrOptim.nl(par=p0, fn=fn, gr=gr, hin=hin, hin.jac=hin.jac) 
ans1

ans2 <- auglag(par=p0, fn=fn, gr=gr, hin=hin, hin.jac=hin.jac) 
ans2

#######################################################################################################
fn <- function(x) (x[1] + 3*x[2] + x[3])^2 + 4 * (x[1] - x[2])^2

gr <- function(x) {
g <- rep(NA, 3)
g[1] <- 2*(x[1] + 3*x[2] + x[3]) + 8*(x[1] - x[2]) 
g[2] <- 6*(x[1] + 3*x[2] + x[3]) - 8*(x[1] - x[2]) 
g[3] <- 2*(x[1] + 3*x[2] + x[3])
g
}

heq <- function(x) {
h <- rep(NA, 1)
h[1] <- x[1] + x[2] + x[3] - 1
h
}


heq.jac <- function(x) {
j <- matrix(NA, 1, length(x))
j[1, ] <- c(1, 1, 1)
j
}

hin <- function(x) {
h <- rep(NA, 1)
h[1] <- 6*x[2] + 4*x[3] - x[1]^3 - 3
h[2] <- x[1]
h[3] <- x[2]
h[4] <- x[3]
h
}


hin.jac <- function(x) {
j <- matrix(NA, 4, length(x))
j[1, ] <- c(-3*x[1]^2, 6, 4)
j[2, ] <- c(1, 0, 0)
j[3, ] <- c(0, 1, 0)
j[4, ] <- c(0, 0, 1)
j
}

set.seed(12)
p0 <- runif(3)
ans <- constrOptim.nl(par=p0, fn=fn, gr=gr, heq=heq, heq.jac=heq.jac, hin=hin, hin.jac=hin.jac) 

set.seed(12)
p0 <- runif(3)
ans2 <- constrOptim.nl(par=p0, fn=fn, heq=heq, hin=hin) 

p0 <- runif(3)
ans3 <- auglag(par=p0, fn=fn, heq=heq, hin=hin) 
#######################################################################################################
##  4-component univariate Gaussian mixture MLE
# Common variance
# 9 parameters to be estimated
#
dgaussmix <- function (p) {
prop <- p[1:nmix]
mu <- p[(nmix+1):(2*nmix)]
sigma <- p[(2*nmix+1)]
sapply(y, function(y)sum(prop*dnorm(y,mean=mu,sd=sqrt(sigma))))	
}

##
rgaussmix <- function (n, prop, mu, sigma) {
nmix <- length(mu)
imix <- sample(1:nmix, size=n, prob=prop, rep=TRUE)
y <- rnorm(n, mean = mu[imix], sd = sqrt(sigma))
return(y)
}


gaussmix.mloglik <- function(p){
sum(log(dgaussmix(p)))
}

gaussmix.grad <- function(p){
g <- rep(NA, length(p))
f <- dgaussmix(p)
pj <- p[1:nmix]
mu <- p[(nmix+1): (2*nmix)]
sigma <- p[2*nmix + 1]
phi <- outer(y, mu, function(y, mu) dnorm(y,mean=mu,sd=sqrt(sigma)))
g[1:nmix] <- - colSums(phi/f)  
phi2 <- outer(y, mu, function(y, mu) (y - mu)/sigma)
fimuj <- t(t(phi * phi2) * pj)
g[(nmix+1): (2*nmix)] <- - colSums(fimuj/f)  
phi3 <- outer(y, mu, function(y, mu) (y - mu)^2/sigma)
fisig <- apply(t(t(phi * ( 1 - phi3) ) * pj), 1, sum)
g[2*nmix+1] <- sum(fisig / f) / (2 * sigma)
-g
}

gmix <- function(p){
par <- rep(NA, length(p)+1)
par[1:3] <- p[1:3]
par[5:(length(p)+1)] <- p[4:length(p)]
par[4] <- 1 - sum(p[1:3])
sum(log(dgaussmix(par)))
}


heq <- function(x) {
x[1] + x[2] + x[3] + x[4] - 1
}


heq.jac <- function(x) {
j <- matrix(NA, 1, length(x))
j[1, ] <- c(1, 1, 1, 1, rep(0, 5))
j
}

hin <- function(x) {
h <- rep(NA, 9)
# all proportions between 0 and 1 
h[1] <- x[1]  
h[2] <- x[2] 
h[3] <- x[3] 
h[4] <- x[4]
h[5] <- 1 - x[1] 
h[6] <- 1 - x[2] 
h[7] <- 1 - x[3] 
h[8] <- 1 - x[4]
# variance greater than 0
h[9] <- x[9] 
h
}


hin.jac <- function(x) {
j <- matrix(NA, 9, length(x))
j[1, ] <- c(1, rep(0,8))
j[2, ] <- c(0, 1, rep(0,7))
j[3, ] <- c(0,0,1, rep(0,6))
j[4, ] <- c(0,0,0,1, rep(0,5))
j[5, ] <- c(-1, rep(0,8))
j[6, ] <- c(0, -1, rep(0,7))
j[7, ] <- c(0,0,-1, rep(0,6))
j[8, ] <- c(0,0,0,-1, rep(0,5))
j[9, ] <- c(rep(0,8), 1)
j
}


p <- c(0.15,0.55,0.2,0.1)
nmix <- length(p)
mu <- c(0,3,6,9)
sigma <- 2
npts <- 1000
y <- rgaussmix(npts, p, mu, sigma)

ymean <- mean(y)
ysd <- sd(y)
p0 <- rep(1/nmix, nmix)
ymean0 <- ymean + ysd * runif(nmix, -1.2, 1.2)
ysd0 <- ysd
par0 <- c(p0,ymean0, ysd0^2)
ans1 <- constrOptim.nl(par=par0, fn=gaussmix.mloglik, gr=gaussmix.grad, heq=heq, heq.jac=heq.jac, hin=hin, hin.jac=hin.jac, control.outer=list(itmax=20), control.optim=list(fnscale=-1))

grad(x=ans1$par[-4], func=gmix)

ans2 <- auglag(par=par0, fn=gaussmix.mloglik, heq=heq, heq.jac=heq.jac, hin=hin, control.outer=list(itmax=20), control.optim=list(fnscale=-1))

####################################################################################################################
# Minimize the following function
fn <- function(x, mat)  {
- prod(c(mat %*% x))
}

hin <- function(par, ...) {
# All constraints are defined such that hin[i] > 0 for all i.
h <- rep(NA, 6)
h[1] <- par[1]
h[2] <- 1 - par[1] 
h[3] <- par[2]
h[4] <- 1 - par[2] 
h[5] <- par[3]
h[6] <- 1 - par[3] 
h
}

heq <- function(par, ...) {
# All equalities are defined such that heq[i] = 0 for all i.
heq <- rep(NA, 1)
heq[1] <- 1 - sum(par)
heq
}

Amat <- matrix(rnorm(9), 3, 3)

project <- function(par, lower, upper) {
# A projection function for generating "feasible" starting values
slack <- 1.e-07
if (par[1] <= 0) par[1] <- slack
if (par[1] >= 1) par[1] <- 1 - slack
if (par[2] <= 0) par[2] <- slack
if (par[3] <= 0) par[3] <- slack
par <- par/sum(par)
par
} 


p0 <- project(runif(3))  # a randomly generated feasible starting value

#p0 <- c(0.1, 0.4, 0.6)
ans <- constrOptim.nl(par=p0, fn=fn, hin=hin, heq=heq, mat=Amat) 
ans
ans2 <- auglag(par=p0, fn=fn, hin=hin, heq=heq, mat=Amat) 
ans2
####################################################################################################################



