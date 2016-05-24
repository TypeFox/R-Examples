# ex1.11.R
set.seed(123)
N <- 100   # number of end points of the grid including T
T <- 1 # length of the interval [0,T] in time units
Delta <- T/N # time increment
W <- numeric(N+1) # initialization of the vector W
t <- seq(0,T, length=N+1)
for(i in 2:(N+1))
        W[i] <- W[i-1] + rnorm(1) * sqrt(Delta)	
x <- 0
y <- -1
BB <- x + W - t/T * (W[N+1] - y +x)
plot(t,BB,type="l")
abline(h=-1, lty=3)
