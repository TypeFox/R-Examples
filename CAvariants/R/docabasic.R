docabasic <-
function (xo,mi,mj) 
{
n<-sum(xo)
 x <- xo/n
#J <- dim(xo)[2]
rsums <- as.vector(rowSums(x))
csums <- as.vector(colSums(x))
dj<-diag(csums)
di<-diag(rsums)
#uni <- matrix(1, 1, J)
drm1 <- diag( 1/( rsums + (rsums==0) ) * (1-(rsums==0)) )
dcm1 <- diag( 1/( csums + (csums==0) ) * (1-(csums==0)) )
drmh <- sqrt(drm1)
dcmh <- sqrt(dcm1)
#ratio <- ( x%*%dcm1 - rsums%*%uni  )*sqrt(n)
R <- drm1 %*% x
C <- dcm1 %*% t(x)
   # rmax <- min(dim(xo)) - 1
    Bpoly <- emerson.poly(mj, csums)$B
 #   Bpoly2 <- sqrt(dj) %*% Bpoly
    Apoly <- emerson.poly(mi, rsums)$B
    #Apoly2 <- (di) %*% Apoly
 #   Z <- t(Apoly)  %*% (ratio)  %*%dj%*% (Bpoly) #useful to check coordinates
Z <- t(Apoly) %*% x  %*% (Bpoly)*sqrt(n) #useful to check coordinates
    ZZ <- Z^2
    pi <- (Apoly) %*% Z %*% t(Bpoly)
    ZtZ <- Z%*%t(Z)
    tZZ <-t(Z)%*%Z 
mu<-diag(ZtZ) 
mu2<-diag(tZZ)
#r<-rmax
#browser()
doca<- new("cabasicresults",
RX=R,CX=C,Cweights=drmh,Rweights=dcmh,Raxes= Bpoly,
Caxes=Apoly,mu=mu,mu2=diag(tZZ),catype="DOCA",tauDen=0,Z=Z,ZtZ=ZtZ,tZZ=tZZ)
}

