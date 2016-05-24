# ex1.06.R
set.seed(123)
N <- 100   # number of end-points of the grid including T
T <- 1 # length of the interval [0,T] in time units
Delta <- T/N # time increment
W <- numeric(N+1) # initialization of the vector W
t <- seq(0,T, length=N+1)
for(i in 2:(N+1))
        W[i] <- W[i-1] + rnorm(1) * sqrt(Delta)	
plot(t,W, type="l", main="Wiener process" , ylim=c(-1,1))
