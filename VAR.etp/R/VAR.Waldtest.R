VAR.Waldtest <-
function(x,p,restrict,type="const")
{
k <- ncol(x)
VARmat = VAR.Rmat(p,k,restrict,type)
Rmat = VARmat$Rmat; rmat = VARmat$rvec
Cmat = VARmat$Cmat; cmat = VARmat$cvec
var1 <- VAR.est(x,p,type)
b <- var1$coef
sigu <- var1$sigu
z <- var1$zmat
b1 = t( t( as.vector(b)) )
tem1 = Cmat %*% b1 - cmat
mat1 <- solve(tcrossprod(z)) %x% sigu
mat2 <- Cmat %*% mat1 %*% t(Cmat)
mat3 <- t(tem1) %*% solve(mat2) %*% tem1
Fstat = as.numeric(mat3/nrow(Cmat))
pval = 1-pf(Fstat,df1=nrow(Cmat),df2=(nrow(x)-p)-ncol(b))
return(list(Fstat=Fstat,pval=pval))
}
