VARB.Boot <-
function(x,p,nb=200,type="const")
{
#set.seed(12345)
n <- nrow(x); k <- ncol(x)
var1 <- VAR.est(x,p,type)
b <- var1$coef
e <- sqrt( (n-p) / ( (n-p)-ncol(b)))*var1$resid


var2 <- VARB.est(x,p,type)
bb <- var2$coef
eb <- sqrt( (n-p) / ( (n-p)-ncol(b)))*var2$resid


mat <- matrix(0,nrow=k,ncol=ncol(b))
matb <- matrix(0,nrow=k,ncol=ncol(b))

for(i in 1:nb)
{
es <- resamp(e)
xs <- VAR.ys(x,b,p,es,type)
bs <- VAR.est(xs,p,type)$coef
mat <- mat + bs/nb

es <- resamp(eb)
xs <- VARB.ys(x,bb,p,es,type)
bs <- VARB.est(xs,p,type)$coef
matb <- matb + bs/nb

}

bias <- mat - b
bs <- VAR.adjust(b,bias,p,type); colnames(bs) <- VAR.names(x,p,type)
es <- VAR.resid(x,bs,var1$zmat,p); colnames(es) <- rownames(b)

biasb <- matb - bb
bsb <- VAR.adjust(bb,biasb,p,type); colnames(bs) <- VAR.names(x,p,type)
esb <- VAR.resid(x,bs,var1$zmat,p); colnames(es) <- rownames(b)


sigu <- t(es) %*% es / ( (n-p) -k*p -1)

return(list(coef=bs,coefb=bsb,resid=es,residb=esb,sigu=sigu,Bias=bias))
}
