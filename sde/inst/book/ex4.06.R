require(sde)
# ex4.06.R
tau0 <- 0.6
k0 <- ceiling(1000*tau0)
set.seed(123)
X1 <- sde.sim(X0=1, N=2*k0, t0=0, T=tau0, model="CIR", theta=c(6,2,1))
X2 <- sde.sim(X0=X1[2*k0+1], N=2*(1000-k0), t0=tau0, T=1, model="CIR", theta=c(6,2,3))

Y <- ts(c(X1,X2[-1]), start=0, deltat=deltat(X1))
X <- window(Y,deltat=0.01) 
DELTA <- deltat(X)
n <- length(X)

mu <- function(x) 6-2*x
sigma <- function(x) sqrt(x)

cpoint(X,mu,sigma)

# nonparametric estimation of the drift
cpoint(X)
