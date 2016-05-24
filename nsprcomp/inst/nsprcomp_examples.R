library(MASS)
set.seed(1)

# Regular PCA, with the tolerance set to return five PCs
prcomp(Boston, tol = 0.36, scale. = TRUE)

# Sparse PCA with different cardinalities per component. The number of components
# is derived from the length of vector k.
nsprcomp(Boston, k = c(13,7,5,5,5), scale. = TRUE)  

# Non-negative sparse PCA with four components. Note that the principal axes
# naturally have a high degree of orthogonality, because each component
# maximizes the additional variance not already explained.
nspc <- nsprcomp(Boston, k = c(7,5,2,2), nneg = TRUE, scale. = TRUE)

# continue the computation of components from a partial model
nsprcomp(Boston, k = 3, ncomp = 5, nneg = TRUE, scale. = TRUE, partial_model = nspc)

# The reconstruction error for each sample can be influenced using the 
# weighting vector omega. To reconstruct the data, the generalized
# inverse of the pseudo-rotation matrix has to be used, because the constrained 
# principal axes are in general not pairwise orthogonal.
X <- matrix(runif(5*10), 5)
nspc <- nsprcomp(X, omega = c(5,1,1,1,5), ncomp = 2, nneg = TRUE)
X_hat <- predict(nspc)%*%ginv(nspc$rotation) + matrix(1,5,1)%*%nspc$center
rowSums((X - X_hat)^2)
