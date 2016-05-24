# Matrix arithmetic
des <- matrix(data = c(9, 5, 8, 4), nrow = 2, byrow = FALSE); des
des + des; des * 2; des * des; 1 / des; des ^ 2; des %% 3
des + 1:3  # warning with differing lengths

# Matrix multiplication
tree <- matrix(data = 1:12, nrow = 3, byrow = FALSE); tree
t(tree) %*% tree
crossprod(tree); tcrossprod(tree)
identical(t(tree) %*% tree, crossprod(tree))

# Matrix inversion
des2 <- solve(des); des2
all.equal(des2 %*% des, des %*% des2)  # TRUE
library(MASS); ginv(des)

# Other matrix functions
upper.tri(x = des, diag = TRUE)  # Upper triangular part: logical values 
des3 <- des                      # Upper triangular part: values 
des3[lower.tri(x = des, diag = FALSE)] <- 0; des3

det(des)                         # determinant
sum(diag(des))                   # trace
qr(des)$rank                     # rank
library(Matrix); rankMatrix(des)
t(des)                           # transpose 
eigen(des)                       # eigenvalue and eigenvector