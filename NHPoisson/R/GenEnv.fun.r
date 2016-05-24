GenEnv.fun <-function(nsim,lambda, fun.name,fun.args=NULL,clevel=0.95,
cores=1,fixed.seed=NULL)
{
cl<-makeCluster(cores)
clusterExport(cl, objects(, envir = .GlobalEnv))
if (!is.null(fixed.seed)) clusterSetRNGStream(cl=cl,iseed=fixed.seed)
simval <- parSapply(cl, c(1:nsim), FUN=funSim.fun,
lambda, fun.name=fun.name,fun.args=fun.args)
stopCluster(cl)


simval<-simval[complete.cases(simval)]
l1<-length(simval)
simval<-matrix(simval, ncol=1)
valmed<-apply(simval,MARGIN=2, FUN=mean)
valinf<-apply(simval,MARGIN=2, FUN=quantile,p=1-clevel) 
valsup<-apply(simval,MARGIN=2, FUN=quantile,p=clevel)

cat('Lower interval: ', valinf, fill=T)
cat('Mean value: ', valmed, fill=T)
cat('Upper interval: ', valsup, fill=T)
cat('Number of valid simulations: ', l1, fill=T)

obj<- list(valmed=valmed,valinf=valinf,valsup=valsup, lambda=lambda, 
	nsim=nsim, nsimval=l1,fixed.seed=fixed.seed)

return(obj)
}
