mfw_eval <- function(d,x,m){

# Computes the multivariate Fourier Whittle criterion for the 
# estimation of the long-memory parameter at value d. 
## 
## The periodogram of data x is defined as
##
##	Ixx(k) = w(x)*conj(w(x)),
##
##	where
##				  N
##	w(k) =  (2*pi*n)^(-1/2)  sum  x(t)*exp(i*2*pi*(k-1)*t/N), 1 <= k <= N.
##                  		 t=1
##
## 	INPUT	d   Long-range dependence parameters (k*1 vector)
##		x   Data (n*k vector)
##		m   Truncation number in Fourier frequencies
##
##	OUTPUT      Fourier Whittle criterion
##
##					based on the paper of Katsumi Shimotsu, 2007
##					Achard & Gannaz (2014)
##_________________________________________________________________________________

if(is.matrix(x)){
	n <- dim(x)[1]
	k <- dim(x)[2]
}else{
	n <- length(x)
	k <- 1
}
x <- as.matrix(x,dim=c(n,k))
if(length(d)!=k){ 
	stop('Number of dependence parameters and number of components in the process mismatch') 
}

t <- seq(0,(n-1),1) 
lambda <- 2*pi*t/n
wx <- matrix(0,n,k)
for(j in 1:k){
    xx <- x[,j]
    wx[,j] <- (2*pi*n)^(-1/2)*Conj(fft(xx))*exp(1i*lambda)
}
wx <- as.matrix(wx[2:(m+1),])

Iwx <- array(0,dim=c(k,k,m));
for(mm in 1:m){
   Iwx[,,mm]<-t(t(wx[mm,]))%*%Conj(wx[mm,])
}

lambda <- lambda[2:(m+1)]
llambda<-matrix(0,m,k)
for(j in 1:k){
	llambda[,j] <- lambda^d[j]*exp((lambda-pi)*d[j]*1i/2)
}
lw <- llambda*wx

g <- t(lw)%*%Conj(lw)/m 
rg <- Re(g)

r <- log(det(rg)) - 2*sum(d)*mean(log(lambda))

return(r)
}
