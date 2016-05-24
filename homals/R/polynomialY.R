`polynomialY` <-
function(d,y,r,verbose=0) {
  k<-length(d)
  zz<-orthogonalPolynomials(d,1:k,min(r,k-1))
  aa<-crossprod(zz,d*y)
  list(y=zz%*%aa,z=zz,a=aa)
}

