# ex1.04.R
set.seed(123)
n <- 10000
beta <-1
K <- 1

x <- rnorm(n)
y <- sapply(x, function(x) max(0,exp(beta*x)-K))

# the true value
exp(beta^2/2)*pnorm(beta-log(K)/beta)-K*pnorm(-log(K)/beta)

t.test(y[1:100]) # first 100 simulations
t.test(y[1:1000]) # first 1000 simulations
t.test(y) # all simulation results

set.seed(123)
x <- rexp(n,rate=0.5)
h <- function(x) (max(0,1-exp(beta*sqrt(x))) + 
  max(0,1-exp(-beta*sqrt(x))))/sqrt(2*pi*x)
y <- sapply(x, h)

# variance reduction
# CALL = PUT + e^{0.5*beta^2} - K
z <- y +exp(0.5*beta^2) - K

t.test(z[1:100]) # first 100 simulations
t.test(z[1:1000]) # first 1000 simulations
t.test(z) # all simulation results
