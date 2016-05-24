###############################################################################
## Simulation study comparing Tukey's biweight with the rmx estimator
###############################################################################

## fixed variables
n <- 11
M <- 1e5
eps.lower <- 0
eps.upper <- 0.05
seed <- 123 ## due to 1e5 replications the influence of the seed is neglectable

eps <- 0.01
contD <- Norm(0, 9)
(res1 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Td(df = 3)
(res2 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Cauchy()
(res3 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(3, 1)
(res4 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(10, 1)
(res5 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(1.51)
(res6 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(1000)
(res7 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))

eps <- 0.02
contD <- Norm(0, 9)
(res11 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Td(df = 3)
(res12 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Cauchy()
(res13 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(3, 1)
(res14 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(10, 1)
(res15 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(1.51)
(res16 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(1000)
(res17 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))

eps <- 0.04
contD <- Norm(0, 9)
(res21 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Td(df = 3)
(res22 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Cauchy()
(res23 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(3, 1)
(res24 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(10, 1)
(res25 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(1.51)
(res26 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(1000)
(res27 <- AffySimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
