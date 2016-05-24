## This script tests 
## - Covariates Coeff Sizes
## - Covariates Coeff Values
## - Mesh Order 2
## written by: Luca Giussani

library(fdaPDE)

order = 1
mesh = create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                      segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5,1)),
                      order = order)
FEMbasis = create.FEM.basis(mesh)

locations = rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0))
observations = c(1,2,1,2,1)
covariates1 = cbind(c(5,2,1,4,5))
covariates2 = cbind(c(5,2,1,4,5), c(5,4,3,2,1))

lambda1 = c(10)
lambda2 = c(10,20)

# SUMMARY
#
# - 1 covariate:
#               |   loc not     |   loc         |
#               |   on nodes    |   on nodes    |
#   ---------------------------------------------
#   1 lambda    |       -       |       -       |
#   ---------------------------------------------
#   2 lambdas   |       -       |       -       |
#   ---------------------------------------------
#
# - 2 covariates:
#               |   loc not     |   loc         |
#               |   on nodes    |   on nodes    |
#   ---------------------------------------------
#   1 lambda    |       -      |      -         |
#   ---------------------------------------------
#   2 lambdas   |       -      |      -         |
#   ---------------------------------------------

output =
    smooth.FEM.basis(observations = observations,
                     locations=locations,
                     FEMbasis = FEMbasis,
                     lambda = lambda1,
                     covariates=covariates1,
                     GCV = TRUE,
                     CPP_CODE = TRUE);
print(output$beta)


output =
    smooth.FEM.basis(observations = observations,
                     locations = locations,
                     FEMbasis = FEMbasis,
                     lambda = lambda1,
                     covariates=covariates2,
                     GCV = TRUE,
                     CPP_CODE = TRUE);
print(output$beta)

output =
    smooth.FEM.basis(observations = observations,
                     locations = locations,
                     FEMbasis = FEMbasis,
                     lambda = lambda2,
                     covariates=covariates1,
                     GCV = TRUE,
                     CPP_CODE = TRUE);
print(output$beta)

output =
    smooth.FEM.basis(observations = observations,
                     locations = locations,
                     FEMbasis = FEMbasis,
                     lambda = lambda2,
                     covariates=covariates2,
                     GCV = TRUE,
                     CPP_CODE = TRUE);
print(output$beta)

output =
    smooth.FEM.basis(observations = observations,
                     FEMbasis = FEMbasis,
                     lambda = lambda1,
                     covariates=covariates1,
                     GCV = TRUE,
                     CPP_CODE = TRUE);
print(output$beta)

output =
    smooth.FEM.basis(observations = observations,
                     FEMbasis = FEMbasis,
                     lambda = lambda1,
                     covariates=covariates2,
                     GCV = TRUE,
                     CPP_CODE = TRUE);
print(output$beta)

output =
    smooth.FEM.basis(observations = observations,
                     FEMbasis = FEMbasis,
                     lambda = lambda2,
                     covariates=covariates1,
                     GCV = TRUE,
                     CPP_CODE = TRUE);
print(output$beta)

output =
    smooth.FEM.basis(observations = observations,
                     FEMbasis = FEMbasis,
                     lambda = lambda2,
                     covariates=covariates2,
                     GCV = TRUE,
                     CPP_CODE = TRUE);
print(output$beta)