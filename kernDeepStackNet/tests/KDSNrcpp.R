library(kernDeepStackNet)
######################
# Check crossprodRcpp

set.seed(100)
A <- matrix(rnorm(100), nrow=20, ncol=5)
AtArcpp <- crossprodRcpp(A)[[1]]
AtA <- base::crossprod(A)
stopifnot(all.equal(AtA, AtArcpp))

##########################
# Check getEigenValuesRcpp

all.equal(rev(getEigenValuesRcpp(AtArcpp)), eigen(AtArcpp)$values)
