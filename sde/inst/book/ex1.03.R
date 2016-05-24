# ex1.03.R
set.seed(123)
n <- 10000
beta <-1
K <- 1

x <- rexp(n,rate=0.5)
h <- function(x) (max(0,1-exp(beta*sqrt(x))) + 
  max(0,1-exp(-beta*sqrt(x))))/sqrt(2*pi*x)
y <- sapply(x, h)

# the true value
K*pnorm(log(K)/beta)-exp(beta^2/2)*pnorm(log(K)/beta-beta) 

t.test(y[1:100]) # first 100 simulations
t.test(y[1:1000]) # first 1000 simulations
t.test(y) # all simulation results