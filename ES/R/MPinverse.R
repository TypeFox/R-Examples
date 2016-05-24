MPinverse <-
function(Matrix)
{
EPS <- 100*(.Machine$double.eps)
U <- eigen(Matrix)$vectors
D <- eigen(Matrix)$values
Dmatrix <- diag(length(D))
Dmatrix[which(Dmatrix>0)] <- 1/D
Dmatrix[which(abs(Dmatrix)>(1/EPS))] <- 0
A <- U%*%Dmatrix%*%t(U)
return(A)
}
