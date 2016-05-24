require(sde)
# ex4.05.R
tau0 <- 0.6
k0 <- ceiling(1000*tau0)
set.seed(123)
X1 <- sde.sim(X0=1, N=2*k0, t0=0, T=tau0, model="CIR", theta=c(6,2,1))
X2 <- sde.sim(X0=X1[2*k0+1], N=2*(1000-k0), t0=tau0, T=1, model="CIR", theta=c(6,2,3))

Y <- ts(c(X1,X2[-1]), start=0, deltat=deltat(X1))
X <- window(Y,deltat=0.01) 
DELTA <- deltat(X)
n <- length(X)

# ex4.05.R (cont)
mu <- function(x) 6-2*x
sigma <- function(x) sqrt(x)
Z <- (diff(X) - mu(X[1:(n-1)])*DELTA)/(sqrt(DELTA)*sigma(X[1:(n-1)]))

tau <- seq(0,1, length=length(Z))
k <- ceiling(n*tau)

Sn <- cumsum(Z^2)
S <- sum(Z^2)
D <- abs((2:n)/n - Sn/S)

k0 <- which(D==max(D))
tau[k0]
sqrt(Sn[k0]/k0)
sqrt((S-Sn[k0])/(n-k0))


# ex4.05.R (cont)
par(mar=c(3,3,1,1))
par(mfrow=c(2,1))
plot(X)
abline(v=tau0,col="red",lty=3)
plot(tau,D,type="l")
abline(v=tau[k0],col="blue")
abline(v=tau0,col="red",lty=3)




