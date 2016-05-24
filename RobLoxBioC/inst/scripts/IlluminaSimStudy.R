###############################################################################
## Simulation study comparing Illumina's default method with the rmx estimator
###############################################################################

library(RobLoxBioC)

## fixed variables
n <- 30
M <- 1e5
eps.lower <- 0
eps.upper <- 0.05
seed <- 123 ## due to 1e5 replications the influence of the seed is neglectable

eps <- 0.01
contD <- Norm(0, 9)
(res1 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Td(df = 3)
(res2 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Cauchy()
(res3 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(3, 1)
(res4 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(10, 1)
(res5 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(3)
(res6 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(1000)
(res7 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))

eps <- 0.02
contD <- Norm(0, 9)
(res11 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Td(df = 3)
(res12 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Cauchy()
(res13 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(3, 1)
(res14 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(10, 1)
(res15 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(3)
(res16 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(1000)
(res17 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))

eps <- 0.04
contD <- Norm(0, 9)
(res21 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Td(df = 3)
(res22 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Cauchy()
(res23 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(3, 1)
(res24 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Norm(10, 1)
(res25 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(3)
(res26 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
contD <- Dirac(1000)
(res27 <- IlluminaSimStudy(n = n, M = M, eps = eps, seed = seed, 
                     eps.lower = eps.lower, eps.upper = eps.upper, 
                     contD = contD))
