enumerateBVS = function(data,forced=NULL,cov=NULL,a1=0,rare=FALSE,mult.regions=FALSE,regions=NULL,hap=FALSE,inform=FALSE){
	p = dim(data)[2]-1
	if(rare==FALSE){
		which.ind = 1:p} 
	if(rare==TRUE & mult.regions==FALSE){
	    which.ind = 1}
	if(rare==TRUE & mult.regions==TRUE){
		which.ind = 1:length(unique(regions))}
	
	#### create all possible models
	model.type=0:1
	all.models <- lapply(vector("list", p), function(v) { model.type } )
	all.models <- expand.grid(all.models)
	
	##Get results for all models
	results <- apply(all.models,1,fitBVS,data=data,forced=forced,cov=cov,a1=a1,rare=rare,mult.regions=mult.regions,
	                 regions=regions,hap=hap,inform=inform)
	coef = t(results[which.ind,])
	if(rare==TRUE & mult.regions==FALSE){
		coef = results[which.ind,]
	}
	fitness = results[(length(which.ind)+1),]
	logPrM = results[length(which.ind)+2,]
	which = all.models
	alpha = rep(a1,dim(all.models)[1])
	which = apply(which,1,paste,collapse="")
	results = list(fitness,logPrM,which,coef,alpha)
	names(results) = c("fitness","logPrM","which","coef","alpha")
	return(results)
	}