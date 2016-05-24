corr.nn4pbo<-function(lam, p, PO.cor=NULL){

x=rnorm(100000,0,1) ##underlying ordinal variable
xcdf=pnorm(x) ##obtain cdf of x

o=rep(0, length(x))
for(i in 1:length(x)){
for(j in 1:length(p)){
if(xcdf[i]>p[j]){o[i]=j}
}}

cxo=cor(x,o) 
cor.xp=PO.cor/cxo

y=rnorm(100000,0,1)
ypois=qpois(pnorm(y),lam)

corrected=cor.xp/cor(ypois,y)
return(corrected)
}