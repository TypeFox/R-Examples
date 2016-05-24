set.seed(123)

#set the number of components K
#set the dimensionality p 
#set the sample size to simulate from
K <- 3
p <- 2
n <- 1000


#sets the parameters to simulate data from 
tau <- c(0.25, 0.3, 0.45)
Mu <- matrix(c(0, 4, 4, 2, 4, 10), 3)
la <- matrix(c(0.2, 0.5, 1, 0.5, 0.5, 0.7), 3)
S <- array(NA, dim = c(p, p, K))
S[,,1] <- matrix(c(0.4, 0, 0, 0.4), 2)
S[,,2] <- matrix(c(4.5, -0.9, -0.9, 2.7), 2)
S[,,3] <- matrix(c(2, -1, -1, 2), 2)


A <- Manly.sim(n, la, tau, Mu, S)
print(A)