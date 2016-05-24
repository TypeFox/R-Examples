library(MASS)
set.seed(1)

# Regular PCA, with tolerance set to return five PCs
prcomp(Boston, tol = 0.35, scale. = TRUE)

# Sparse cumulative PCA with five components and 20 non-zero loadings. 
# The orthonormality penalty is set to a value which avoids co-linear principal 
# axes. Note that the non-zero loadings are not distributed uniformly over 
# the components.
nscumcomp(Boston, ncomp = 5, k = 20, gamma = 1e3, scale. = TRUE)  

# Non-negative sparse cumulative PCA
nscumcomp(Boston, ncomp = 5, nneg = TRUE, k = 20, gamma = 1e3, scale. = TRUE)  
