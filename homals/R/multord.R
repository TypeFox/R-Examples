`multOrd` <-
function(d,y,r,itermax=100,eps=1e-6,verbose=0){
z<-orthogonalPolynomials(d,1:length(d),r)
a<-crossprod(z,d*y)
z1<-cbind(z[,1]); a1<-rbind(a[1,])
z2<-cbind(z[,2:r]); a2<-rbind(a[2:r,])
iter<-1; sold<-Inf
repeat{
	ytilde<-y-z2%*%a2
	z1<-tcrossprod(ytilde,a1)/sum(a1^2)
	z1<-twoDirections(z1,d); z1<-z1/sqrt(sum(d*z1^2))
	a1<-crossprod(z1,d*ytilde)
	ytilde<-y-z1%*%a1
  	qq<-La.svd(sqrt(d)*ytilde,r-1,r-1)
  	z2<-(1/sqrt(d))*qq$u
	z<-weightedGramSchmidt(cbind(z1,z2),d)$pol     
	a<-crossprod(z,d*y); a2<-rbind(a[2:r,])
  	snew<-sum(d*(y-z%*%a)^2)
	if (verbose > 2) print(c(round(sold,6),round(snew,6)))
	if ((iter == itermax) || ((sold - snew) < eps) || (snew < eps)) break()
	    iter<-iter+1; sold<-snew
	}
  list(yhat=z%*%a,z=z,a=a)
}

