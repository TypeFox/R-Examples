CL.nsteps <- function (model, fun, npoints, nsteps, lower, upper, parinit = NULL, kmcontrol = NULL, control=NULL) {
    
    n <- nrow(model@X)

    if (length(kmcontrol$penalty)==0) kmcontrol$penalty <- model@penalty
    if (length(kmcontrol$optim.method)==0) kmcontrol$optim.method <- model@optim.method
    if (length(kmcontrol$parinit)==0) kmcontrol$parinit <- model@parinit
    if (length(kmcontrol$control)==0) kmcontrol$control <- model@control
    kmcontrol$control$trace=FALSE

    for (i in 1:nsteps) {
		L <- min(model@y)
	  	res.CL <- max_qEI.CL(model, npoints=npoints, L=L, lower=lower, upper=upper, parinit=parinit, control=control)
	  	
	  	model@X <- rbind(model@X, res.CL$par)     
      model@y <- rbind(model@y, as.matrix(apply(res.CL$par, 1, fun), npoints))
	  	
	  	model <- km(formula=model@trend.formula, design=model@X, response=model@y, 
	  		    covtype=model@covariance@name, lower=model@lower, upper=model@upper, 
                            nugget=NULL, penalty=kmcontrol$penalty, optim.method=kmcontrol$optim.method, 
		    	    parinit=kmcontrol$parinit, control=kmcontrol$control, gr=model@gr)
	  	            
	}

    return(list(par = model@X[(n + 1):(n + nsteps*npoints),, drop = FALSE], 
                value = model@y[(n + 1):(n + nsteps*npoints),, drop = FALSE], 
                npoints=npoints, nsteps=nsteps, lastmodel = model))
}



