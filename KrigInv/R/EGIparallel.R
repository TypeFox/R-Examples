EGIparallel <- function(T, model, method=NULL, method.param=NULL, 
						fun, iter, batchsize = 1, lower, upper, new.noise.var=0, 
						optimcontrol=NULL, kmcontrol=NULL,integcontrol=NULL,...) {
  
	n <- nrow(model@X); d <- model@d
	
	if (is.null(kmcontrol$penalty)) kmcontrol$penalty <- model@penalty
	if (length(model@penalty==0)) kmcontrol$penalty <- NULL 
	if (is.null(kmcontrol$optim.method)) kmcontrol$optim.method <- model@optim.method 
	if (is.null(kmcontrol$parinit)) kmcontrol$parinit <- model@parinit
	if (is.null(kmcontrol$control)) kmcontrol$control <- model@control
	if (is.null(kmcontrol$CovReEstimate)) kmcontrol$CovReEstimate <- model@param.estim
	if (is.null(method)) method <- "sur"
	if (is.null(optimcontrol$optim.option)) optimcontrol$optim.option <- 2
  
	for (i in 1:iter) {
		if (method == "sur"){
			
			integration.param <- integration_design(integcontrol,d,lower,upper,model,T)
			
			oEGI <- max_sur_parallel(lower=lower, upper=upper,optimcontrol=optimcontrol,T=T,model=model,
								integration.param=integration.param,batchsize=batchsize,
								new.noise.var=new.noise.var)
		}
    else if(method == "timse" || method == "imse"){
      integration.param <- integration_design(integcontrol,d,lower,upper,model,T)
      
      oEGI <- max_timse_parallel(lower=lower, upper=upper,optimcontrol=optimcontrol,T=T,model=model,
  							integration.param=integration.param,batchsize=batchsize,
								new.noise.var=new.noise.var,epsilon=method.param,imse=(method == "imse"))
    }
		else if(method == "vorob"){
		  integration.param <- integration_design(integcontrol,d,lower,upper,model,T)
		  
		  oEGI <- max_vorob_parallel(lower=lower, upper=upper,optimcontrol=optimcontrol,T=T,model=model,
		                           integration.param=integration.param,batchsize=batchsize,
		                           new.noise.var=new.noise.var)
      
		}
       		
		print("New points"); print(oEGI$par)
		X.new <- oEGI$par; y.new <- rep(0,times=nrow(X.new))
		for (i in 1:nrow(X.new)) y.new[i] <- fun((oEGI$par)[i,],...)
		
		model <- update_km(model=model,NewX=X.new,NewY=y.new,CovReEstimate=kmcontrol$CovReEstimate,new.noise.var=rep(new.noise.var,times=batchsize),kmcontrol=kmcontrol)		
	}
	
	return(list(
				par=model@X[(n+1):(n+iter*batchsize),, drop=FALSE], 
				value=model@y[(n+1):(n+iter*batchsize),, drop=FALSE], 
				npoints=batchsize, 
				nsteps=iter*batchsize, 
				lastmodel=model,
				lastvalue=oEGI$value,
				allvalues=oEGI$allvalues
				)
			)
}
