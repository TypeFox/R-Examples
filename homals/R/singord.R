`singOrd` <-
function(d,y,r,itermax=100,eps=1e-6,verbose=0){
z<-cbind(1:length(d)); z<-z-sum(d*z)/sum(d); z<-z/wNorm(d,z)
a<-crossprod(z,d*y)
iter<-1; sold<-Inf
repeat{
	z<-tcrossprod(y,a)/sum(a^2)
	z<-twoDirections(z,d); z<-cbind(z/wNorm(d,z))
	a<-crossprod(z,d*y)
  	snew<-sum(d*(y-z%*%a)^2)
	if (verbose > 2) print(c(round(sold,6),round(snew,6)))
	if ((iter == itermax) || ((sold - snew) < eps) || (snew < eps)) break()
		iter<-iter+1; sold<-snew
	}
list(yhat=z%*%a,z=z,a=a)
}

