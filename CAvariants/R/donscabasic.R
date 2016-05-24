donscabasic <-
function (xo,mi,mj) 
{
#  rmax <- min(dim(xo)) - 1
x <- xo/sum(xo)
rsums <- as.vector(rowSums(x))
csums <- as.vector(colSums(x))
drm1 <- diag( 1/( rsums + (rsums==0) ) * (1-(rsums==0)) )
dcm1 <- diag( 1/( csums + (csums==0) ) * (1-(csums==0)) )
drmh<-diag(rep(1,nrow(x))) #change the metric in NSCA
dcmh <- sqrt(dcm1)
dj <- diag(csums)
di <- diag(rsums)
tauden<-1 - sum(rsums^2)
#Apoly <- emerson.poly(mi, rsums)$BT  #with trivial
Apoly <- emerson.poly(mi, rsums)$B  #without trivial
Apoly2 <- sqrt(di) %*% Apoly
Bpoly <- emerson.poly(mj, csums)$B
Bpoly2 <- sqrt(dj) %*% Bpoly
#pcc <- 1/sqrt(tauden)*(drmh %*% ( x - rsums %*% t(csums) ) %*% dcm1)
pcc <- (drmh %*% ( x - rsums %*% t(csums) ) %*% dcm1)
 Z <- t(Apoly2) %*% pcc %*% dj %*% Bpoly #no trivial
  pi <- (Apoly2) %*% Z %*% t(Bpoly2) #no trivial
ZtZ<-Z%*%t(Z)   
 mu <- diag(ZtZ)
#tau<-sum(mu)/tauden
tZZ<-t(Z)%*%Z
    mu2<- diag(tZZ)
#r<-rmax
    Cweights <- dj
#browser()
donsca<- new("cabasicresults",
RX=pcc,CX=t(pcc),Rweights=dj,Cweights=diag(rep(1,nrow(xo))),
        Raxes=Bpoly,Caxes=Apoly2,mu=mu,mu2=mu2,catype="DONSCA",tauDen=tauden,Z=Z,ZtZ=ZtZ,tZZ=tZZ)
}

