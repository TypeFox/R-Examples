
par.est <- function(data=data) 
{
  
if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data) < 2)  stop("data must have at least 2 columns")
  
n <- nrow(data)
p <- ncol(data)-1
  
Y <- data[, 1]
X <- data[, -1]
if (!is.matrix(X))  X <- as.matrix(X) 
  
X.X <- t(X)%*%X
X.Y <- t(X)%*%Y
beta.est <- symsolve(X.X, X.Y)
  
return(beta.est)
  
}   
