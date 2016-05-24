pls_eigen <-
function(X,Y,a){

# Eigenvector method for PLS

X <- as.matrix(X)
Y <- as.matrix(Y)
P=eigen(t(X)%*%Y%*%t(Y)%*%X)$vectors[,1:a]
T=X%*%P
Q=eigen(t(Y)%*%X%*%t(X)%*%Y)$vectors[,1:a]
U=Y%*%Q

list(P=P,T=T,Q=Q,U=U)
}

