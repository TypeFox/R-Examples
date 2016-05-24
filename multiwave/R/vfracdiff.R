vfracdiff <- function(x,d){

## Vectorial fractional differencing procedure
##
## INPUT   x	Initial process
##	   d	Vector of long-range dependence parameters
##
## OUTPUT 	diag((1-L)^d) x
##
##				Achard & Gannaz (2014)
##___________________________________________________________

if(is.matrix(x)){
	N <- dim(x)[1]
	k <- dim(x)[2]
}else{
	N <- length(x)
	k <- 1
}
x <- as.matrix(x,dim=c(N,k))
dx <- x

if(length(d)!=k){
	stop('The number of components in x and the number of dependence parameters d must be the same')
}

a <- seq(1,N-1,1)
for(i in 1:k){
	b <- (a+d[i]-1)/a
	b <- c(1,cumprod(b))
	temp <- fft( fft(c(x[,i],rep(0,N))) * fft(c(b,rep(0,N))), inverse=TRUE )/(2*N)
	dx[,i] <- Re(temp[1:N])
}

if(k==1){
dx <- as.vector(dx)
}

return(dx)

}
