gcontrol<-function(data,zeta=1000,kappa=4,tau2=1,epsilon=0.01,ngib=500,burn=50,idum=2348)
{
  nkdata <- dim(data)[1]
  deltot <- rep(0,nkdata)
  kdata <- as.matrix(data)
  x <- a <- rep(0,nkdata)
  z<-.C("gcontrol",kdata=as.double(kdata),nkdata=as.integer(nkdata),
        zeta=as.double(zeta),kappa=as.double(kappa),tau2=as.double(tau2),
        epsilon=as.double(epsilon),ngib=as.integer(ngib),burn=as.integer(burn),
        idumR=as.integer(idum),deltot=as.double(deltot),x=as.double(array(x)),A=as.double(a),PACKAGE="gap")

  nkdata6 <- nkdata/6
  list(deltot=z$deltot[1:nkdata6],x2=z$x[1:nkdata6],A=z$A[1:nkdata6])
}
