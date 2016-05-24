## This script tests 
## - PDE smoothing 
## - 2nd order FEs 
## - C++ code

library(fdaPDE)

data(mesh.2D.simple)
observations = sin(pi*mesh.2D.simple$nodes[,1]) + rnorm(n = nrow(mesh.2D.simple$nodes), sd = 0.1)

FEMbasis = create.FEM.basis(mesh.2D.simple)

lambda = c(10^-2, 10^-1, 0.5, 5, 10)

PDE_parameters_anys = list(K = matrix(c(0.01,0,0,1), nrow = 2), b = c(0,0), c = 0)
FEM_CPP_PDE = smooth.FEM.PDE.basis(observations = observations, 
                                   FEMbasis = FEMbasis, 
                                   lambda = lambda, 
                                   PDE_parameters = PDE_parameters_anys)
print(FEM_CPP_PDE$fit.FEM$coeff)