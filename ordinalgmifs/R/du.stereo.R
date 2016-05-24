du.stereo <-
function(Xb, x, alpha, phi, Ymat, k) {
		phi <- c(1, phi)
     	eta<-matrix(0,ncol=k-1,nrow=dim(x)[1])
     	for (i in 1:(k-1)) {
            eta[,i]<- exp(alpha[i] + phi[i]*Xb)
        }
		denom <- 1 + .Call("do_matrix_sum_rows", eta)
      	numer<-matrix(0,ncol=k-1,nrow=dim(x)[1])
      	for (i in 1:(k-1)){
      		numer[,i]<-eta[,i]*phi[i]
      	}
      	numer2<- -(.Call("do_matrix_sum_rows",numer))*Ymat[,k]/denom     #contribution to log-like for class K
      	numer1<-matrix(0,ncol=k-1,nrow=dim(x)[1])       # contribution to log-like for classes 1 to K-1
      	for (i in 1:(k-1)){
      		numer1[,i]<-Ymat[,i]*(phi[i]-apply(numer,1,sum)/denom)
      	}
		numer1f <- .Call("do_matrix_sum_rows", numer1) 
      	dll<-x*(numer1f + numer2)  
      	u<-colSums(-dll, na.rm=TRUE)   # u is a 2p*1 vector of -dlogL/dbeta
		update.value<-min(u)
		update.j<-which.min(u)
		list(update.j=update.j, update.value=update.value)
}
