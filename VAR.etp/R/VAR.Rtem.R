VAR.Rtem <-
function(x,p,sigu,Rmat,rmat,type="const",method="gls")
{

n <- nrow(x); k <- ncol(x)
y <- t(x[(p+1):n,])
y1 <- t( t( as.vector(y)) )

z0 <- t(x[n:1,])
z <- matrix(1,nrow=k*p,ncol=1)
for( i in n:(p+1) )
{
    m <- t( t( as.vector(matrix(z0[,(i-p+1):i])) ) )
    z <- cbind(z,m)
}

if(type=="none") z <- z[,2:(n-p+1),drop=FALSE]
if(type=="const") z <- rbind(z[,2:(n-p+1)],matrix(1,nrow=1,ncol=(n-p)))
if(type=="const+trend") z <- rbind(z[,2:(n-p+1)],matrix(1,nrow=1,ncol=(n-p)),matrix((p+1):n,ncol=(n-p)))

z1 <- y1 - (t(z) %x% diag(k)) %*% rmat

if (method=="ols") {
mat1 <- tcrossprod(z) %x% diag(k)
mat2 <- t(Rmat) %*% ( z %x% diag(k) ) %*% z1
br1 <- solve( t(Rmat) %*% mat1 %*% Rmat) %*% mat2
b <- matrix(Rmat %*% br1 +rmat,nrow=k); e <- y - b%*%z; rownames(b) <- colnames(x); colnames(b) <- VAR.names(x,p,type)}
if (method=="gls") {
mat1 <- tcrossprod(z) %x% solve(sigu)
mat2 <- t(Rmat) %*% ( z %x% solve(sigu) ) %*% z1
br2 <- solve( t(Rmat) %*% mat1 %*% Rmat) %*% mat2
b <- matrix(Rmat %*% br2 +rmat,nrow=k); e <- y - b%*%z; rownames(b) <- colnames(x); colnames(b) <- VAR.names(x,p,type)}

sigu <- tcrossprod(e) / (n-p )
tem1 = tcrossprod(z) %x% solve(sigu)
tem2 = Rmat %*% solve( t(Rmat) %*% (tem1) %*%  Rmat) %*% t(Rmat);
tem3 = sqrt(diag(tem2));
tem4 = matrix(tem3,nrow=k,ncol=ncol(b))
tmat = b/tem4; tmat[is.nan(tmat)] = NA 
 
return(list(coef=b,resid=t(e),sigu=sigu,zmat=z,tratio=tmat))
}
