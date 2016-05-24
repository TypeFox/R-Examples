plot.fclust <-
function (x, v1v2, colclus, umin, ucex, pca,...)  
{
fclust.obj=x
if ((missing(fclust.obj)) || (!inherits(fclust.obj, "fclust")))
stop("An object of class fclust must be given")
X=fclust.obj$X
Xca=fclust.obj$Xca
k=fclust.obj$k
n=nrow(X)
p=ncol(X)
U=fclust.obj$U
H=fclust.obj$H
stand=fclust.obj$stand
is.wholenumber<-function(x, tol=.Machine$double.eps^0.5) {abs(x-round(x))<tol}
if (missing(pca))
{
  pca=FALSE
}
if (!is.logical(pca)) 
{
  pca=FALSE
  cat("pca is not logical: the default value pca=FALSE will be used ",fill=TRUE)
}
if (missing(v1v2))
{
v1v2=c(1,2)
}
else
{
if (is.null(v1v2))
{
v1v2=c(1,2)
}
}
if (length(v1v2)!=2)
{
v1v2=c(1,2)
cat("v1v2 must be a vector with two elements specifying the numbers of the variables to be plotted: the default value v1v2=c(1,2) will be used ",fill=TRUE)
}
if (!is.numeric(v1v2))
{
v1v2=c(1,2)
cat("v1v2 must contain integers: the default value v1v2=c(1,2) will be used ",fill=TRUE)
}  
else
{
if (all(is.wholenumber(v1v2))==FALSE)
{
v1v2=c(1,2)
cat("v1v2 must contain integers: the default value v1v2=c(1,2) will be used ",fill=TRUE)
}
if (pca==TRUE)
{
if ((min(v1v2)<1) || (max(v1v2)>min(n,p)))
{
v1v2=c(1,2)
if (n>=p)
cat("v1v2 must contain integers in {1, 2, ..., p}: the default value v1v2=c(1,2) will be used ",fill=TRUE)
else
cat("v1v2 must contain integers in {1, 2, ..., n}: the default value v1v2=c(1,2) will be used ",fill=TRUE)
}
}
else
{  
if ((min(v1v2)<1) || (max(v1v2)>p))
{
v1v2=c(1,2)
cat("v1v2 must contain integers in {1, 2, ..., p}: the default value v1v2=c(1,2) will be used ",fill=TRUE)
}
}
}
if (missing(umin))
{
  umin=0
}
if (!is.numeric(umin)) 
{
  umin=0
  cat("umin is not numeric: the default value umin=0 will be used ",fill=TRUE)
}
if ((umin<0) || (umin>1))  
{
  umin=0
  cat("umin must be in the interval [0,1]: the value umin=0 will be used ",fill=TRUE)
}
if (missing(ucex))
{
  ucex=FALSE
}
if (!is.logical(ucex)) 
{
  ucex=FALSE
  cat("ucex is not logical: the default value ucex=FALSE will be used ",fill=TRUE)
}
if (missing(colclus))
{
palette(rainbow(k))  
}
else
{
if (is.null(colclus))
{
palette(rainbow(k))  
}
else
{
palette(colclus)  
if (length(colclus)<k)
{
  cat("The length of colclus is lower than the number of clusters k: points belonging to different clusters will have same colors ",fill=TRUE)
}
}  
}
info.U=cl.memb.t(U,umin)
colunit=c()
for (i in 1:n)
{
  if (info.U[i,1]==0)
  {
  colunit[i]="black"
  }
else
  {
  colunit[i]=palette()[info.U[i,1]]
  }
}
if (ucex==TRUE)
{
dimunit=info.U[,2]+0.5
dimunit[which(info.U[i,1]==0)]=0.5
}
else
{
dimunit=1.5
}  
if (pca==FALSE)
{
nx=colnames(X)[v1v2[1]]  
ny=colnames(X)[v1v2[2]]  
if (stand==1)
H=Hraw(X,H)
plot(X[,v1v2[1]],X[,v1v2[2]],xlab=nx,ylab=ny,pch=16,col=colunit,cex=dimunit)
points(H[,v1v2[1]],H[,v1v2[2]],col=1:k,pch=8,cex=2,lwd=2)
}
else
{
nx=paste("Principal Component",v1v2[1])
ny=paste("Principal Component",v1v2[2])
if (stand==0)
{
Xca=scale(X,center=TRUE,scale=TRUE)[,]
H=(H-matrix(1,k,p)%*%diag(apply(X,2,mean),nrow=p))/(matrix(1,k,p)%*%diag(apply(X,2,sd),nrow=p))
}
s=svd(Xca)
sc.unit=Xca%*%s$v
sc.cent=H%*%s$v
fitpc=(s$d[v1v2[1]]+s$d[v1v2[2]])/sum(s$d)*100
plot(sc.unit[,v1v2[1]],sc.unit[,v1v2[2]],xlab=nx,ylab=ny,pch=16,col=colunit,cex=dimunit,sub=paste0("Variability explained by these two components: ",round(fitpc,2), "%"))
points(sc.cent[,v1v2[1]],sc.cent[,v1v2[2]],col=1:k,pch=8,cex=2,lwd=2)
}
}