		# Compute g-inverse of information matrix
invInfMat <-
function(C,N,T){

ei <- eigen(t(T) %*% t(N) %*% C %*% N %*% T)
nn <- length(ei$values)
L <- matrix(0, nrow=nn, ncol=nn)
for(i in 1:(nn)) {
if( Re(ei$values[i]) <1e-6)next
L <- L + (1/Re(ei$values[i]))*Re(ei$vectors[,i])%*%t(Re(ei$vectors[,i]))
}
return(L)
}
