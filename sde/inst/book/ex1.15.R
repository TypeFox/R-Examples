# ex1.15.R
W <- BM()
t <- time(W)
N <- length(t)
x <- 10
theta <- 5
sigma <- 3.5
X <- numeric(N)
X[1] <- x
ito.sum <- c(0, sapply(2:N, function(x) { 
   exp(-theta*(t[x]-t[x-1])) * (W[x]-W[x-1])} ) )
X <- sapply(1:N, function(x)  {X[1]*exp(-theta*t[x]) + 
    sum(ito.sum[1:x])} )
X <- ts(X,start=start(W), deltat=deltat(W))
plot(X,main="Ornstein-Uhlenbeck process")
