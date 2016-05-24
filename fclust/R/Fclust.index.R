Fclust.index <-
function (fclust.obj, index, alpha)
{
if ((missing(fclust.obj)) || (!inherits(fclust.obj,"fclust")))
stop("An object of class fclust must be given")
X=fclust.obj$Xca
U=fclust.obj$U
H=fclust.obj$H
m=fclust.obj$m
lg=length(c("PC","PE","MPC","SIL","SIL.F","XB","ALL"))
if (missing(index))
{
indexN=lg
}
else
{
if (is.null(index))
{
indexN=lg
}
else
{
indexN=match(toupper(index),c("PC","PE","MPC","SIL","SIL.F","XB","ALL"))
}
}
if (any(is.na(indexN)))
{
indexN=lg
cat("(At least one) no match is found: all the indices will be computed ",fill=TRUE)
}
index=indexN
out.index=rep(0,lg-1)
if (any(index==1) || any(index==lg))
out.index[1]=PC(U)
if (any(index==2) || any(index==lg))
out.index[2]=PE(U)
if (any(index==3) || any(index==lg))
out.index[3]=MPC(U)
if (any(index==4) || any(index==lg))
out.index[4]=SIL(X,U)$sil
if (any(index==5) || any(index==lg))
{ 
if (missing(alpha))
{
alpha=1
cat("The default value alpha=1 has been set for computing SIL.F ",fill=TRUE)
}
if (!is.numeric(alpha)) 
{
alpha=1
cat("The weighting coefficient alpha is not numeric: the default value alpha=1 will be used for computing SIL.F ",fill=TRUE)
}
if (alpha<0)  
{
alpha=1
cat("The number of clusters k must be non negative: the value alpha=1 will be used for computing SIL.F ",fill=TRUE)
}
out.index[5]=SIL.F(X,U,alpha)
} 
if (any(index==6) || any(index==lg))
out.index[6]=XB(X,U,H,m)
names(out.index)=c("PC","PE","MPC","SIL","SIL.F","XB")
if (max(index)<lg)
out.index=out.index[index]
return(out.index)
}
