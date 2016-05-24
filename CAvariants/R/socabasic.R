socabasic <-
function (xo,mj) 
{
n<-sum(xo)
    x <- xo/n
rsums <- as.vector(rowSums(x))
csums <- as.vector(colSums(x))
di<-diag(rsums)
dj<-diag(csums)
    Bpoly <- emerson.poly(mj, csums)$B
#  Bpoly <- orthopoly.exe(c(csums))[,-1]

    Bpoly2 <- sqrt(dj) %*% Bpoly
######################################################
drm1 <- diag( 1/( rsums + (rsums==0) ) * (1-(rsums==0)) )
dcm1 <- diag( 1/( csums + (csums==0) ) * (1-(csums==0)) )
drmh <- sqrt(drm1)
dcmh <- sqrt(dcm1)
ratio <- drmh %*% ( x - rsums %*% t(csums) ) %*% dcmh*sqrt(n)
u<-svd(ratio)$u
mu <- svd(ratio)$d
R <- drm1 %*% x
C <- dcm1 %*% t(x)
Z <- t(u) %*% ratio %*% Bpoly2 #useful to check coordinates
ZtZ <- Z%*%t(Z)
tZZ <- t(Z)%*%Z
#rmax <- min(dim(xo)) - 1
#r<-rmax
soca <- new("cabasicresults",
          RX=R,CX=C,Rweights=dcmh,Cweights=drmh,
         Raxes= Bpoly,Caxes=u,mu=mu,mu2=diag(tZZ),tauDen=0,catype="SOCA",Z=Z,ZtZ=ZtZ,tZZ=tZZ)

}
