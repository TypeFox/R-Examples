sonscabasic <-
function (xo,mj) 
{
#rmax <- min(dim(xo)) - 1
x <- xo/sum(xo)
rsums <- as.matrix(rowSums(x))
csums <- as.vector(colSums(x))
tauden <- 1 - sum(rsums^2)
drm1 <- diag( 1/( rsums + (rsums==0) ) * (1-(rsums==0)) )
dcm1 <- diag( 1/( csums + (csums==0) ) * (1-(csums==0)) )
drmh<-diag(rep(1,nrow(x))) #change the metric in NSCA
dcmh <- sqrt(dcm1)
dj <- diag(csums)
di <- diag(rsums)
uni <- matrix(1, 1, ncol(x))
uni1 <- rep(1, nrow(x))
Bpoly <- emerson.poly(mj, csums)$B
Bpoly2 <- sqrt(dj) %*% Bpoly
#pcc <- 1/sqrt(tauden)* ( x%*%dcm1 - rsums %*% (uni) ) 
pcc <- ( x%*%dcm1 - rsums %*% (uni) ) 

u<-svd(pcc%*%sqrt(dj))$u

#Z <- t(u)  %*%pcc %*% sqrt(dj)%*%Bpoly2
Z <- t(u)  %*%pcc %*% dj%*%Bpoly

ZtZ<-Z%*%t(Z)
tZZ<-t(Z)%*%Z
mu2<- diag(tZZ) #only the sum gives me the total inertia
mu<- diag(ZtZ) #these are coincident with each eigenvalue (mu^2)
#tau<-sum(mu)
#r<-rmax
#browser()
sonsca <- new("cabasicresults",
          RX=pcc,CX=t(pcc),Rweights=dj,Cweights=diag(uni1),
          Raxes=Bpoly,Caxes=u,mu=mu,mu2=mu2,catype="SONSCA",tauDen=tauden,Z=Z,ZtZ=ZtZ,tZZ=tZZ)

#sonsca <- list(RX=pcc,CX=t(pcc),Rweights=dj,Cweights=diag(uni1),
#          Raxes=Bpoly,Caxes=u,mu=mu,mu2=mu2,tauDen=tauden,catype="SONSCA",Z=Z,ZtZ=ZtZ,tZZ=tZZ)

#class(sonsca)<-"cabasicresults"
#getClass("cabasicresults")
}
