VAR.Cest <-
function(x,p,type="const",Cmat, cvec)
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
mat1 = solve(tcrossprod(z));
b=tcrossprod(y,z) %*% mat1; rownames(b) <- colnames(x); colnames(b) <- VAR.names(x,p,type)
e <- y - b%*%z; 
sigu <- tcrossprod(e) / ( n-p)

mat2 = mat1 %x% sigu;
mat3 = t(Cmat) %*% solve(Cmat %*% mat2 %*% t(Cmat))%*%(cvec-Cmat %*% t(t(as.vector(b))));
br = t(t(as.vector(b))) + mat2 %*% mat3;
br = matrix(br,nrow=k,ncol=ncol(b));
e <- y - br%*%z; 
sigu <- tcrossprod(e) / ( n-p)

return(list(coef=br,resid=t(e),sigu=sigu))
}
