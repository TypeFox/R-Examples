EGI <- function(T, model, method=NULL, method.param=NULL, 
            fun, iter, lower, upper, new.noise.var=0,
            optimcontrol=NULL, kmcontrol=NULL,integcontrol=NULL,...) {

	n <- nrow(model@X); d <- model@d
	
	if (is.null(kmcontrol$penalty)) kmcontrol$penalty <- model@penalty
	if (length(model@penalty==0)) kmcontrol$penalty <- NULL 
	if (is.null(kmcontrol$optim.method)) kmcontrol$optim.method <- model@optim.method 
	if (is.null(kmcontrol$parinit)) kmcontrol$parinit <- model@parinit
	if (is.null(kmcontrol$control)) kmcontrol$control <- model@control
	if (is.null(kmcontrol$CovReEstimate)) kmcontrol$CovReEstimate <- model@param.estim
	if (is.null(method)) method <- "ranjan"
	
	for (i in 1:iter) {
	
    if (method == "timse" || method == "imse"){
			
			integration.param <- integration_design(integcontrol,d,lower,upper,model,T)			
			oEGI <- max_timse(lower=lower, upper=upper,optimcontrol=optimcontrol,T=T,model=model,
					integration.param=integration.param,
					new.noise.var=new.noise.var,epsilon=method.param,imse=(method == "imse"))
			
		}
		else if (method == "sur" || method == "jn"){
			
			integration.param <- integration_design(integcontrol,d,lower,upper,model,T)
			real.volume.variance <- (method == "jn")
			
			oEGI <- max_sur(lower=lower, upper=upper,optimcontrol=optimcontrol,T=T,model=model,
						integration.param=integration.param,
						new.noise.var=new.noise.var,real.volume.variance=real.volume.variance)
		}
    else{
			oEGI <- max_infill_criterion(lower=lower, upper=upper,T=T, method=method, 
								method.param=method.param,model=model,optimcontrol=optimcontrol)
		}
    
		print("New point"); print(oEGI$par)
		X.new <- oEGI$par; y.new <- fun(oEGI$par,...)
		X.new <- as.numeric(X.new); X.new <- matrix(X.new,nrow=1,ncol=d)
		
		model <- update_km(model=model,NewX=X.new,NewY=y.new,CovReEstimate=kmcontrol$CovReEstimate,new.noise.var=new.noise.var,kmcontrol=kmcontrol) #call to update_km
		
	}
	return(list(
				par=model@X[(n+1):(n+iter),, drop=FALSE], 
				value=model@y[(n+1):(n+iter),, drop=FALSE], 
				npoints=1, 
				nsteps=iter, 
				lastmodel=model,
				lastvalue=oEGI$value,
				allvalues=oEGI$allvalues,
				variance.volume=oEGI$variance.volume
				)
			)
}
