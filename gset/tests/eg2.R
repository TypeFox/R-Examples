# Example 2:  calculation of nminmax
library(gset)

L <- -0.2
U <- 0.2
theta <- 0
sigma <- 0.4  
alpha <- 0.05
beta  <- 0.2
K <- 3
r <- 1

# the sample size per group with a traditional nonsequential design
n.fix <- nfix(r, L,U,theta,sigma,alpha,beta)

# nminmax with nonbinding futility
bound1  <- nminmax(L, U, theta, sigma, n.fix$n1, n.fix$n2, 1:K/K, alpha, beta,n.sim=100)
figureEF(bound1, K)

