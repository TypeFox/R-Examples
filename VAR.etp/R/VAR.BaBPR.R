VAR.BaBPR <-
function(x,p,h,nboot=500,nb=200,type="const",alpha=0.95)
{
k <- ncol(x)
VAR <- VARB.Boot(x,p,nb,type)
b1 <- VAR$coef;  e1 <- VAR$resid
b2 <- VAR$coefb; e2 <- VAR$residb
BIAS <- VAR$Bias
f <- VAR.Fore(x,b1,p,h)

fstat <- array(0,c(h,k,nboot))
for(i in 1:nboot){
es <- resamp(e2)
ys <- VARB.ys(x,b2,p,es,type)
vars <- VAR.est(ys,p)
bs <- vars$coef
bss <- VAR.adjust(bs,BIAS,p,type)
fs <- VAR.ForeB(x,bss,p,h,e1)
fstat[ , ,i] <- as.matrix(fs)
}
PIs <- VAR.PI(fstat,alpha,f)$PI
return(list(Intervals=PIs,Forecast=f,Prob=alpha))
}
