VCV <-
function (Xca, U, H, which)
{
if (missing(Xca))
stop("The data set must be given")
if (is.null(Xca))
stop("The data set Xca is empty")
n=nrow(Xca)
Xca=as.matrix(Xca)
if (any(is.na(Xca)))
stop("The data set Xca must not contain NA values")
if (!is.numeric(Xca))
stop("The data set Xca is not a numeric data.frame or matrix")
if (is.null(rownames(Xca)))
rn=paste("Obj",1:n,sep=" ")
else
rn=rownames(Xca)
if (missing(U))
stop("The membership degree matrix U must be given")
if (is.null(U))
stop("The membership degree matrix U is empty")
U=as.matrix(U)
if (any(is.na(U)))
stop("The membership degree matrix U must not contain NA values")
if (!is.numeric(U))
stop("The membership degree matrix U is not numeric")
if (missing(H))
  stop("The prototype matrix H must be given")
if (is.null(H))
  stop("The prototype matrix H is empty")
H=as.matrix(H)
if (any(is.na(H)))
  stop("The prototype matrix H must not contain NA values")
if (!is.numeric(H)) 
  stop("The prototype matrix H is not numeric")
if (nrow(U)!=nrow(Xca)) 
  stop("The numbers of rows of U and Xca must be the same")
if (nrow(H)!=ncol(U)) 
  stop("The number of rows of H and the one of columns of U must be the same")
if (ncol(H)!=ncol(Xca)) 
  stop("The numbers of columns of H and Xca must be the same")
if (ncol(U)==1) 
  stop("There is only k=1 cluster: VCV is not given")
is.wholenumber<-function(x, tol=.Machine$double.eps^0.5) {abs(x-round(x))<tol}
if (missing(which))
{
  which=c(1,2)
}
else
  {
  if (is.null(which))
  {
    which=c(1,2)
  }
}
if (length(which)>2)
{
  which=c(1,2)
  cat("which must contain no more than two elements specifying the requested plots: the default value which=c(1,2) will be used ",fill=TRUE)
}
if (!is.numeric(which))
{
  which=c(1,2)
  cat("which must contain 1 and/or 2: the default value which=c(1,2) will be used ",fill=TRUE)
}  
else
{
  if (all(is.wholenumber(which))==FALSE)
  {
    which=c(1,2)
    cat("which must contain 1 and/or 2: the default value which=c(1,2) will be used ",fill=TRUE)
  }  
  if ((all(is.wholenumber(which))==TRUE) && ((any(which<1)) || (any(which>2))))
  {
    which=c(1,2)
    cat("which must contain 1 and/or 2: the default value which=c(1,2) will be used ",fill=TRUE)
  }
}  
show <- rep(FALSE, 2)
show[which] <- TRUE
# VAT image
if (show[1])
{
if (length(which)==2)  
  devAskNewPage(ask=TRUE)
VAT(Xca)
}
# VCV image
if (show[2])
{
if (length(which)==2)  
  devAskNewPage(ask=TRUE)
k=nrow(H)
D=as.matrix(dist(H))
setI1=c(1,rep(0,k-1))
Pk=c(1,rep(0,k-1))
for (c in 2:k)
{
 Dr=matrix(D[setI1>0,],nrow=sum(setI1))
 Dr[,setI1==1]=max(Dr)
 wm=(which(Dr==min(Dr), arr.ind = TRUE))[2]
 setI1[wm]=1
 Pk[c]=wm
}
U=U[,Pk]
Ucm=cl.memb(U)
U=U[order(Ucm[,1],-Ucm[,2]),]
Xca=Xca[order(Ucm[,1],-Ucm[,2]),]
#U=U[order(Ucm[,1],Ucm[,2]),]
Rstar=matrix(0,nrow=n,ncol=n)
for (i1 in 1:n-1)
{
for (i2 in (i1+1):n)
{
di1i2c=rep(0,k)
for (c in 1:k)
{  
di1i2c[c]=sqrt(sum((Xca[i1,]-H[c,])^2))+sqrt(sum((Xca[i2,]-H[c,])^2))
}  
Rstar[i1,i2]=min(di1i2c)  
Rstar[i2,i1]=Rstar[i1,i2]
}
}
image(Rstar[ncol(Rstar):1,],col=grey(seq(0, 1, length = 256)),xlab="",ylab="",xaxt="n",yaxt="n",main="VCV")
}
if (length(which)==2)  
  devAskNewPage(ask=FALSE)
}