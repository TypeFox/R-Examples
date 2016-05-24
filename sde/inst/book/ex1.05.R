# ex1.05.R
set.seed(123)
n <- 10000
beta <-1
K <- 1
x <- rnorm(n)
y1 <- sapply(x, function(x) max(0,K-exp(beta*x)))
y2 <- sapply(-x, function(x) max(0,K-exp(beta*x)))

y <- (y1+y2)/2
# the true value
K*pnorm(log(K)/beta)-exp(beta^2/2)*pnorm(log(K)/beta-beta) 

t.test(y[1:100]) # first 100 simulations
t.test(y[1:1000]) # first 1000 simulations
t.test(y) # all simulation results
