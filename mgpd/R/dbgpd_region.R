dbgpd_region <-
function(x,y,z,quant=c(.95,0.9,.75),...){
prec              = 1000
sq                = (x[2]-x[1])*(y[2]-y[1])
one               = sum(z)*sq
level             = seq(quantile(z,0.4,na.rm=TRUE),quantile(z,0.9995,na.rm=TRUE),length.out = prec)
numint            = rep(0,prec)
levelf            = function(l) sum(z[z>l])*sq
numint            = sapply(level,levelf)/one
nq                = rep(0,length(quant))
for(i in 1:length(quant)) nq[i] = level[numint<=quant[i]][1]
out               = list(one,quant,nq,x,y,z)
names(out)        = c("one","prob","q","x","y","z")
out
}
