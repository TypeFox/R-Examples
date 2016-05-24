working.fitvariog1 <-
function(gammadat,models=c('expn','mten','sphn','gaun'))
{
 ## gammadat : x : distance, gamma empirical semivariance
 ## models: variogram models to fit
 ## estimated coefficient of a better valid model
 pii <- with(gammadat,firstpeak(x,gamma))
 range0 <- with(gammadat,x[pii])
 ii <- which(gammadat$n>0)
 gammadat <- gammadat[ii,]
 out <- working.compvariogmodels1(gammadat,models,range0)
 if(max(out$parms)==0){ret=F}
 else{ret=T}
 list(coef=c(out$parms[2],range0,out$parms[1]),model=out$model,ret=ret)
}
