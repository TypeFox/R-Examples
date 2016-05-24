VIFCR <-
function (fclust.obj, which)  
{
if ((missing(fclust.obj)) || (!inherits(fclust.obj, "fclust")))
stop("An object of class fclust must be given")
Xca=fclust.obj$Xca
n=nrow(Xca)
p=ncol(Xca)
U=fclust.obj$U
H=fclust.obj$H
k=fclust.obj$k
is.wholenumber<-function(x, tol=.Machine$double.eps^0.5) {abs(x-round(x))<tol}
if (missing(which))
{
  which=1:3
}
else
{
  if (is.null(which))
  {
    which=1:3
  }
}
if (length(which)>3)
{
  which=1:3
  cat("which must contain no more than three elements specifying the requested plots: the default value which=1:3 will be used ",fill=TRUE)
}
if (!is.numeric(which))
{
  which=1:3
  cat("which must contain integers in {1,2,3}: the default value which=1:3 will be used ",fill=TRUE)
}  
else
{
  if (all(is.wholenumber(which))==FALSE)
  {
    which=1:3
    cat("which must contain integers in {1,2,3}: the default value which=1:3 will be used ",fill=TRUE)
  }
  if ((all(is.wholenumber(which))==TRUE) && ((any(which<1)) || (any(which>3))))
{
    which=1:3
    cat("which must contain integers in {1,2,3}: the default value which=1:3 will be used ",fill=TRUE)
  }
}
show <- rep(FALSE, 3)
show[which] <- TRUE
# PLOT 1
if (show[1])
{
if (length(which)>1)  
    devAskNewPage(ask=TRUE)
vU=as.vector(U)
rangeU=seq(0,1,0.1)
cardU=c()
sumU=c()
for (lev in 1:9){
  cardU[lev]=length(which((vU>=rangeU[lev]) & (vU<rangeU[lev+1])))
  sumU[lev]=sum(vU[which((vU>=rangeU[lev]) & (vU<rangeU[lev+1]))])
}
cardU[10]=length(which((vU>=rangeU[10]) & (vU<=rangeU[11])))
sumU[10]=sum(vU[which((vU>=rangeU[10]) & (vU<=rangeU[11]))])
hU=hist(vU,breaks=seq(0,1,0.1),plot=FALSE)
hU$counts=((((k*(k-2))/(k-1))*sumU)+((k/(k-1))*cardU))/n
plot(hU,main="Cluster Balance",xlab="",ylab="",yaxt='n')
}
# PLOT 2
if (show[2])
{
if (length(which)>1)  
  devAskNewPage(ask=TRUE)
x=c()
y=c()
for (i in 1:n)
{
  x[i]=max(U[i,])
  y[i]=max(U[i,][-(which.max(U[i,]))])
}
plot(x,y,main="Cluster Max Memb. Degrees",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")  
}
# PLOT 3
if (show[3])
{
if (length(which)>1)  
  devAskNewPage(ask=TRUE)
if (k==2) par(mfrow=c(1,2))
if ((k==3) || (k==4)) par(mfrow=c(2,2))
if ((k==5) || (k==6)) par(mfrow=c(2,3))
if ((k==7) || (k==8)) par(mfrow=c(2,4))
if (k==9) par(mfrow=c(3,3))
if (k>9){ 
  di=k%/%4
  if (k%%4>0) di=di+1
  par(mfrow=c(di,4))
}
D=matrix(0,nrow=n,ncol=k)
for (i in 1:n) 
{
  for (c in 1:k) 
  {
    D[i,c]=sum((Xca[i,]-H[c,])^2)
  }
}
for (c in 1:k)
plot(D[,c],U[,c],main=paste("Cluster",c),xlab="",ylab="",xaxt='n',yaxt='n')  
}
par(mfrow=c(1,1))
if (length(which)>1)  
  devAskNewPage(ask=FALSE)
}