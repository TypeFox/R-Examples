SIL <-
function (Xca, U)
{
if (missing(Xca))
stop("The data set must be given")
if (is.null(Xca))
stop("The data set in Xca is empty")
n=nrow(Xca)
Xca=as.matrix(Xca)
if (any(is.na(Xca)))
stop("The data set in Xca must not contain NA values")
if (!is.numeric(Xca)) 
stop("The data set in Xca is not a numeric data.frame or matrix")
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
stop("There is only k=1 cluster: the silhouette index is not computed")
if (nrow(U)!=nrow(Xca)) 
stop("The numbers of rows of U and Xca must be the same")
memb=vector(length(n),mode="numeric")
for (i in 1:n)
{
memb[i]=which.max(U[i,])
}
k=ncol(U)
count=c()
for(c in 1:k)
{
count[c]=length(which(memb==c))
c=c+1
}
D=matrix(0,nrow=n,ncol=n)
for (i in 1:(n-1))
{
for (i2 in (i+1):n)
{
D[i,i2]=sum((Xca[i,]-Xca[i2,])^2)
D[i2,i]=D[i,i2]
}
}
a=rep(0,n)
b=a
sil.obj=a
B=matrix(0,nrow=n,ncol=k)
for (i in 1:n)
{
for (c in 1:k)
{
for(i2 in 1:n)
{
if (memb[i2]==c) 
 B[i,c]=B[i,c]+D[i,i2]
}
}
}
for (i in 1:n)
{
for (c in 1:k)
{
if (memb[i]==c)
{
if (count[c]!=1)
{
B[i,c]=B[i,c]/(count[c]-1)
a[i]=B[i,c]	
B[i,c]=max(B[i,])+1	
}
}
else 
B[i,c]=B[i,c]/count[c]
}
if (count[memb[i]]!=1)
{
b[i]=min(B[i,])
sil.obj[i]=(b[i]-a[i])/max(a[i],b[i])
}
else
sil.obj[i]=0
}
names(sil.obj)=rn
sil=mean(sil.obj)
out=list()
out$sil.obj=sil.obj
out$sil=sil
return(out)
}
