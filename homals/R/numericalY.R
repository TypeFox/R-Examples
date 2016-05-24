`numericalY` <-
function(d,y,r,verbose=0) {
  z0<-orthogonalPolynomials(d,1:length(d),1)
  a0<-as.vector(crossprod(z0*d,y))
  if (r == 1)
  	return(list(y=z0%o%a0,z=cbind(z0),a=rbind(a0)))
  else {
  	yy<-y-z0%o%a0
  	qq<-La.svd(sqrt(d)*yy,r,r)
  	zz<-cbind(z0,array((1/sqrt(d))*qq$u[,1:(r-1)],c(dim(y)[1],r-1)))
  	aa<-rbind(a0,array((qq$d[1:r]*qq$vt)[1:(r-1),],c(r-1,dim(y)[2])))
  	return(list(y=zz%*%aa,z=zz,a=aa))}
}

