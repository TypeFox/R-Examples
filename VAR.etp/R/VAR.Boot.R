VAR.Boot <-
function(x,p,nb=200,type="const")
{
#set.seed(12345)
n <- nrow(x); k <- ncol(x)
var1 <- VAR.est(x,p,type)
b <- var1$coef
e <- sqrt( (n-p) / ( (n-p)-ncol(b)))*var1$resid

mat <- matrix(0,nrow=k,ncol=ncol(b))
for(i in 1:nb)
{
es <- resamp(e)
xs <- VAR.ys(x,b,p,es,type)
bs <- VAR.est(xs,p,type)$coef
mat <- mat + bs/nb
}

bias <- mat - b
bs <- VAR.adjust(b,bias,p,type); colnames(bs) <- VAR.names(x,p,type)
es <- VAR.resid(x,bs,var1$zmat,p); colnames(es) <- rownames(b)
sigu <- t(es) %*% es / ( (n-p) -ncol(b))

return(list(coef=bs,resid=es,sigu=sigu,Bias=bias))
}
