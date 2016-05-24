cabasic <-
function(X) { 
n<-sum(X)
X <- X/n
#rmax <- min(dim(X))-1
rsums <- as.vector(rowSums(X))
csums <- as.vector(colSums(X))
di<-diag(rsums)
dj<-diag(csums)
drm1 <- diag( 1/( rsums + (rsums==0) ) * (1-(rsums==0)) )
dcm1 <- diag( 1/( csums + (csums==0) ) * (1-(csums==0)) )
drmh <- sqrt(drm1)
dcmh <- sqrt(dcm1)
#ratio <- (drmh %*% ( X - rsums %*% t(csums) ) %*% dcmh)*n
ratio <- drmh %*% ( X - rsums %*% t(csums) ) %*% dcmh
ratio2<-drm1%*%X%*%dcm1
Y <- svd(ratio)
mu <- Y$d
#r <- sum(mu>1e-15)
#r<-rmax
RX <- drm1 %*% X
CX <- dcm1 %*% t(X)
#if (r < rmax) { 
 # mu[ (r+1):rmax ] <- 0
 # Raxes[ , (r+1):rmax ] <- 0
  #Caxes[ , (r+1):rmax ] <- 0
#}

#ca<-list(RX=RX,CX=CX,Rweights=dcmh,Cweights=drmh,
#          Raxes=Y$v,Caxes=Y$u,r=r,mu=mu,mu2=0,catype="CA",tau=0,tauDen=0,Z=ratio2,ZtZ=RX,tZZ=RX)

#setClass("cabasicresults",
#representation(
#  RX="matrix", CX="matrix", Rweights="matrix", Cweights="matrix",
#  Raxes="matrix", Caxes="matrix", r="numeric", mu="numeric",mu2="numeric",catype="character",
#tau="numeric",tauDen="numeric",Z="matrix",ZtZ="matrix",tZZ="matrix"),S3methods=FALSE )

cabasic<-new("cabasicresults", RX=RX,CX=CX,Rweights=dcmh,Cweights=drmh,
          Raxes=Y$v,Caxes=Y$u,mu=mu,mu2=0,catype="CA",tauDen=0,Z=ratio2,ZtZ=RX,tZZ=RX)

#class(ca)<-"cabasicresults"
cabasic
}

