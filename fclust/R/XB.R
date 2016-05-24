XB <-
function (Xca, U, H, m)
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
stop("There is only k=1 cluster: the XB index is not computed")
k=ncol(U)
if (missing(m))
{
m=2
}
if (!is.numeric(m)) 
{
m=2
cat("The parameter of fuzziness m is not numeric: the default value m=2 will be used ",fill=TRUE)
}
if (m<=1) 
{
m=2
cat("The parameter of fuzziness m must be >1: the default value m=2 will be used ",fill=TRUE)
}
D=matrix(0,nrow=n,ncol=k)
for (i in 1:n) 
{
for (c in 1:k) 
{
D[i,c]=sum((Xca[i,]-H[c,])^2)
}
}
distH=10^10*sum(H^2)
for (c1 in 1:(k-1)) 
{
for (c2 in (c1+1):k) 
{
if (sum((H[c1,]-H[c2,])^2)<distH)
distH=sum((H[c1,]-H[c2,])^2)
}	
}
xie.beni=sum((U^m)*D)/(n*distH)
return(xie.beni)
}
