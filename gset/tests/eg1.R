# example 1:  equivalence boundary only
library(gset)
  
r <- 1
theta <- 0
L <- -0.2
U <- 0.2
sigma <- 0.4  
alpha <- 0.05
beta  <- 0.2
K <- 3
      
# the sample size per group with a traditional nonsequential design
n.fix <- nfix(r, L,U,theta,sigma,alpha,beta)
      
      
# default 
# there are two ways to generate the boundary plots
# 1. specify plot=TRUE (default) in "binding()"
equivonly(L, U, sigma, n.fix$n1, n.fix$n2, 1:K/K, alpha, plot=FALSE, n.sim=100)             
      
# 2. specify plot=FALSE in "binding()" and apply the "figureE()" command 
bound  <- equivonly(L, U,  sigma, n.fix$n1, n.fix$n2, 1:K/K, alpha, plot=FALSE, n.sim=100)  
figureE(bound, K)
