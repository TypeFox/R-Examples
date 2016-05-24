VAR.est <-
function(x,p,type="const")
{
n <- nrow(x); k <- ncol(x)
y <- t(x[(p+1):n,])
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

b <- tcrossprod(y,z) %*% solve(tcrossprod(z)); rownames(b) <- colnames(x); colnames(b) <- VAR.names(x,p,type)
e <- y - b%*%z; 
sigu <- tcrossprod(e) / ( (n-p)-ncol(b))
zz <- tcrossprod(z) /(n-p)
tem1 = (n-p)^(-1) * solve(zz)  %x% sigu
tem2 = sqrt(diag(tem1))
tem3 = matrix(tem2,nrow=k,ncol=ncol(b))
tmat = b/tem3
return(list(coef=b,resid=t(e),sigu=sigu,zzmat=zz,zmat=z,tratio=tmat))
}
