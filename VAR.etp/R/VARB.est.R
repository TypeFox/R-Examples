VARB.est <-
function(x,p,type="const")
{
n <- nrow(x); k <- ncol(x)
y <- t(x[1:(n-p),])
z0 <- t(x)

z <- matrix(1,nrow=k*p,ncol=1)
for( i in 1:(n-p) )
{
    m <- t( t( as.vector(matrix(z0[,(i+1):(p+i)])) ) )
    z <- cbind(z,m)
}

if(type=="none") z <- z[,2:(n-p+1),drop=FALSE]
if(type=="const") z <- rbind(z[,2:(n-p+1)],matrix(1,nrow=1,ncol=(n-p)))
if(type=="const+trend") z <- rbind(z[,2:(n-p+1)],matrix(1,nrow=1,ncol=(n-p)),matrix((p+1):n,ncol=(n-p)))

b <- tcrossprod(y,z) %*% solve(tcrossprod(z)); rownames(b) <- colnames(x); colnames(b) <- VAR.names(x,p,type)
e <- y - b%*%z; rownames(e) <- rownames(b); colnames(e) <- NULL; 
sigu <- tcrossprod(e) / ( (n-p)-ncol(b))
zz <- tcrossprod(z) /(n-p)
return(list(coef=b,resid=t(e),sigu=sigu,zzmat=zz,zmat=z))
}
