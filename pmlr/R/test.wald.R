# p covariates; J categories
# Wald test for H_0: C * vecB = 0, df = J (ith covariate)
# Ainv = inverse of the Fisher information matrix
# vecB = pJ column vector (i.e., [b_{1,1}, b{1,2}, ..., b{1,J}, b_{2,1}, b_{2,2}, ..., b_{2,J}, ..., b_{p,1}, b{p,2}, ..., b{p,J}] 
# C =  J x pJ constrast matrix

test.wald <- function(vecB, Ainv, C) {
   t(C %*% vecB) %*% solve(C %*% Ainv %*% t(C)) %*% (C %*% vecB)
}
