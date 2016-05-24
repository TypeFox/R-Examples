VAR.Wald <-
function(x,p,restrict,type="const",bootstrap=0,nb=500)
{
k <- ncol(x)
M=VAR.Waldtest(x,p,restrict,type)
Fstat = M$Fstat; pval = M$pval

M1= VAR.Rest(x,p,restrict,type,method="gls")
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
Fmat[i,]=VAR.Waldtest(xs,p,restrict,type)$Fstat
}
Bpval= mean(as.numeric((Fmat > as.numeric(Fstat))))}

return(list(Fstat=Fstat,pval=pval,Boot.pval=Bpval))
}
