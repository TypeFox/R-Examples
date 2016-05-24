`nominalY` <-
function(d,y,r,verbose=0) {
  qq<-La.svd(sqrt(d)*y,r,r)
  zz<-(1/sqrt(d))*qq$u
  aa<-qq$d[1:r]*qq$vt
  list(y=zz%*%aa,z=zz,a=aa)
}

