## This script tests 
## - PDE smoothing with SV coefficients 
## - 1st order FEs
## - C++ code

library(fdaPDE)

data(mesh.2D.rectangular)
observations = sin(0.2*pi*mesh.2D.rectangular$nodes[,1]) + rnorm(n = nrow(mesh.2D.rectangular$nodes), sd = 0.1)

FEMbasis = create.FEM.basis(mesh.2D.rectangular)

lambda = c(10^-2)

K_func<-function(points)
{
  mat<-c(0.01,0,0,1)
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i] = 0.5*mat %*% t(points[i,1]^2)
  output
}
b_func<-function(points)
{
  output = array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 0
  output
}

c_func<-function(points)
{
  rep(c(0), nrow(points))
}

u_func<-function(points)
{
  rep(c(0), nrow(points))
}
# Space-varying smoothing
PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
FEM_CPP_PDE_SV = smooth.FEM.PDE.sv.basis(observations = observations, 
                                      FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters)
print(FEM_CPP_PDE_SV$fit.FEM$coeff)
