getBeta = function(object, burnin = 0, thin = 1){
	if(all(!is.na(object$draws$Beta))){	
	xx = object$draws$Beta
	nn = dim(xx)[[1]]
	dd=seq((burnin+1), nn, thin)
	if(length(dim(xx)) == 3){ #for random effect model
		return(xx[dd,,]) }
	if(length(dim(xx)) == 2){  ##for fixed effect model
		return(xx[dd,]) }
	}else(return(NULL))
}

getIntercept = function(object, burnin = 0, thin = 1){
	xx = object$draws$Intercept
	if(class(xx) == 'matrix'){  ##for fixed effect model
		nn = dim(xx)[[1]]
		dd = seq((burnin+1), nn, thin)
		return(xx[dd,]) 
	}else{			#for random effect model
		nn = length(xx)
		dd = seq((burnin+1), nn, thin)
		return(xx[dd])
	}
}

getAlpha = function(object, burnin = 0, thin = 1){
	if(all(!is.na(object$draws$Alpha))){
	xx = object$draws$Alpha
	nn = length(xx)
	dd = seq((burnin+1), nn, thin)
	return(xx[dd]) 
	}else(return(NULL))
}

getLS = function(object, burnin = 0, thin = 1){
	xx = object$draws$ZZ
	nn = length(xx)
	dd = seq((burnin+1),nn,thin)
	kk = length(xx[[1]])
	lp = list()
	for(i in 1:kk){
		lp.sub = array(0,dim=c(dim(xx[[1]][[i]])[1],dim(xx[[1]][[i]])[2],length(dd)))
		for(j in 1:length(dd)){
			ind = dd[j]
			lp.sub[,,j] = xx[[ind]][[i]]
	}
	lp[[i]] = lp.sub
	}
	return(lp)
#	return(sapply(dd,function(w) lapply(1:kk, function(y) xx[[w]][[y]]))
}	

getLikelihood = function(object, burnin = 0, thin = 1){
	xx = object$draws$likelihood 
	nn = length(xx)
	dd = seq(burnin, nn, thin)
	return(xx[dd]) 
}


