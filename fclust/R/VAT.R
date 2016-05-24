VAT <-
function (Xca)
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
# VAT image
setI1=rep(0,n)
P=c()
R=as.matrix(dist(Xca))
wm=which(R==max(R), arr.ind = TRUE)[2]
setI1[wm]=1
P[1]=wm
for (i in 2:n)
{
Rr=matrix(R[setI1>0,],nrow=sum(setI1))
Rr[,setI1==1]=max(Rr)
wm=(which(Rr == min(Rr), arr.ind = TRUE))[2]
setI1[wm]=1
P[i]=wm
}
Rstar=R[P,rev(P)]
#x=seq(0,length=nrow(Rstar)+1)
#image(x,x,Rstar,col = grey(seq(0, 1, length = 256)),xaxp=c(1,nrow(Rstar),1),yaxp=c(1,nrow(Rstar),1))
image(Rstar,col=grey(seq(0, 1, length = 256)),xlab="",ylab="",xaxt="n",yaxt="n",main="VAT")
}