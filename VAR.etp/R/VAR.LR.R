VAR.LR <-
function(x,p,restrict0,restrict1,type="const",bootstrap=0,nb=500)
{
k <- ncol(x)
M=VAR.LRtest(x,p,restrict0,restrict1,type)
LR = M$LRstat; pval = M$LRpval

M1= VAR.Rest(x,p,restrict0,type,method="gls")
b1=M1$coef; e=M1$resid
nop= sum(as.numeric(b1 !=0))
e1=sqrt( nrow(e)/(nrow(e)-nop))* e

if(bootstrap==0) Bpval=NULL
if(bootstrap>0){
Fmat = matrix(NA,nrow=nb); 
for (i in 1:nb){
if(bootstrap==1) es = resamp(e1)
if(bootstrap==2) es = Mammen(nrow(e1)) * e1
xs=VAR.ys(x,b1,p,es,type)
Fmat[i,]=VAR.LRtest(xs,p,restrict0,restrict1,type)$LRstat
}
Bpval= mean(as.numeric((Fmat > as.numeric(LR))))}

return(list(LRstat=LR,pval=pval,Boot.pval=Bpval))
}
