VCV2 <-
function (Xca, U, which)
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
if (ncol(U)==1)
stop("There is only k=1 cluster: VCV2 is not given")
if (nrow(U)!=nrow(Xca))
stop("The numbers of rows of U and Xca must be the same")
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
# VCV2 image
if (show[2])
{
if (length(which)==2)  
  devAskNewPage(ask=TRUE)
setI1=rep(0,n)
P=c()
U=matrix(1,ncol=n,nrow=n)-(U%*%t(U)/max(U%*%t(U)))
wm=which(U==max(U), arr.ind = TRUE)[2]
setI1[wm]=1
P[1]=wm
for (i in 2:n)
{
Ur=matrix(U[setI1>0,],nrow=sum(setI1))
Ur[,setI1==1]=max(Ur)
wm=(which(Ur == min(Ur), arr.ind = TRUE))[2]
setI1[wm]=1
P[i]=wm
}
Ustar=U[P,rev(P)]
image(Ustar,col=grey(seq(0, 1, length = 256)),xlab="",ylab="",xaxt="n",yaxt="n",main="VCV2")
}
if (length(which)==2)  
  devAskNewPage(ask=FALSE)
}