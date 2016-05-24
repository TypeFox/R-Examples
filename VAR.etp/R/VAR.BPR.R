VAR.BPR <-
function(x,p,h,nboot=500,type="const",alpha=0.95)
{
k <- ncol(x)
var1 <- VAR.est(x,p)
var2 <- VARB.est(x,p,type)
b1 <- var1$coef; e1 <- var1$resid
b2 <- var2$coef; e2 <- var2$resid
f <- VAR.Fore(x,b1,p,h)

fstat <- array(0,c(h,k,nboot))
for(i in 1:nboot){
es <- resamp(e2)
ys <- VARB.ys(x,b2,p,es,type)
vars <- VAR.est(ys,p)
bs <- vars$coef;
fs <- VAR.ForeB(x,bs,p,h,e1)
fstat[ , ,i] <- as.matrix(fs)
}
PIs <- VAR.PI(fstat,alpha,f)$PI
return(list(Intervals=PIs,Forecast=f,Prob=alpha))
}
