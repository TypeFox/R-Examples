X <- model.matrix(model); X             
# hat matrix
X %*% solve(t(X) %*% X) %*% t(X)  
