# ex1.10.R
set.seed(123)
r <- 1
sigma <- 0.5
x <- 10
N <- 100   # number of end points of the grid including T
T <- 1 # length of the interval [0,T] in time units
Delta <- T/N # time increment
W <- numeric(N+1) # initialization of the vector W
t <- seq(0,T, length=N+1)
for(i in 2:(N+1))
        W[i] <- W[i-1] + rnorm(1) * sqrt(Delta)	
S <- x * exp((r-sigma^2/2)*t + sigma*W)
plot(t,S,type="l",main="geometric Brownian motion")
