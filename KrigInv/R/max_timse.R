
max_timse <- function(lower, upper, optimcontrol=NULL, 
		integration.param=NULL,T, model, new.noise.var=0,epsilon=0,imse=FALSE){
	
	if(is.null(integration.param)) integration.param <- integration_design(integcontrol=NULL,d=model@d,
				lower=lower,upper=upper,model=model,T=T)
	if(is.null(epsilon)) epsilon <- 0
	integration.points <- as.matrix(integration.param$integration.points)
	integration.weights <- integration.param$integration.weights
	
	d <- model@d
	if(is.null(optimcontrol$method)) optimcontrol$method<-"genoud"
	
	#precalculates the kriging mean and variance on the integration points
	pred <- predict_nobias_km(object=model,newdata=integration.points,type="UK",se.compute=TRUE)
	intpoints.oldmean<-pred$mean ; intpoints.oldsd<-pred$sd 
	
	#precompute other important data
	precalc.data <- precomputeUpdateData(model,integration.points)
	
	# COMPUTE WEIGHT BEFORE OPTIMIZATION
	if(!imse){
		weight <- 1/sqrt(2*pi*(intpoints.oldsd^2+epsilon^2)) * exp(-0.5*((intpoints.oldmean-T)/sqrt(intpoints.oldsd^2+epsilon^2))^2)
		weight[is.nan(weight)] <- 0
	}else{
		weight <- NULL
	}
	#end of precomputations
	########################################################################################
	#discrete Optimisation
	if(optimcontrol$method=="discrete"){
		if (is.null(optimcontrol$optim.points)){
			n.discrete.points<-d*100
			optimcontrol$optim.points <- t(lower + t(matrix(runif(d*n.discrete.points),ncol=d)) * (upper - lower))
		}
		optim.points <- optimcontrol$optim.points
		optim.points<-data.frame(optim.points)
		all.crit <- seq(1,nrow(optim.points))
		
		for (i in 1:nrow(optim.points)){
			all.crit[i] <- timse_optim(x=t(optim.points[i,]), integration.points=integration.points,integration.weights=integration.weights,
					intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
					precalc.data=precalc.data,
					T=T, model=model, new.noise.var=new.noise.var,weight=weight)
		}
			
		ibest <- which.min(all.crit)
		o <- list(3)
		o$par <- optim.points[ibest,]
		o$value <- min(all.crit)
		o$allvalues<-all.crit
		
		o$value <- as.matrix(o$value)
		colnames(o$par) <- colnames(model@X)
		colnames(o$value) <- colnames(model@y)   
		return(list(par=o$par, value=o$value,allvalues=o$allvalues))
	}
	
	########################################################################################
	#Optimization with Genoud
	if(optimcontrol$method=="genoud"){
		if (is.null(optimcontrol$pop.size))  optimcontrol$pop.size <- 50*d#floor(4 + 3 * log(d))
		if (is.null(optimcontrol$max.generations))  optimcontrol$max.generations <- 10*d#100*d
		if (is.null(optimcontrol$wait.generations))  optimcontrol$wait.generations <- 2#2
		if (is.null(optimcontrol$BFGSburnin)) optimcontrol$BFGSburnin <- 2#10#0
		if (is.null(optimcontrol$parinit))  optimcontrol$parinit <- NULL
		if (is.null(optimcontrol$unif.seed))  optimcontrol$unif.seed <- 1
		if (is.null(optimcontrol$int.seed))  optimcontrol$int.seed <- 1
		if (is.null(optimcontrol$print.level))  optimcontrol$print.level <- 1
    
		#mutations
		if (is.null(optimcontrol$P1)) optimcontrol$P1<-0#50		#copy
		if (is.null(optimcontrol$P2)) optimcontrol$P2<-0#50
		if (is.null(optimcontrol$P3)) optimcontrol$P3<-0#50
		if (is.null(optimcontrol$P4)) optimcontrol$P4<-0#50
		if (is.null(optimcontrol$P5)) optimcontrol$P5<-50
		if (is.null(optimcontrol$P6)) optimcontrol$P6<-50#50
		if (is.null(optimcontrol$P7)) optimcontrol$P7<-50
		if (is.null(optimcontrol$P8)) optimcontrol$P8<-50
		if (is.null(optimcontrol$P9)) optimcontrol$P9<-0
		
		domaine <- cbind(lower, upper)
		
		o <- genoud(fn=timse_optim, nvars=d, max=FALSE, pop.size=optimcontrol$pop.size,
				max.generations=optimcontrol$max.generations,wait.generations=optimcontrol$wait.generations,
				hard.generation.limit=TRUE, starting.values=optimcontrol$parinit, MemoryMatrix=TRUE,
				Domains=domaine, default.domains=10, solution.tolerance=0.000000001,
				boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE, BFGS=TRUE,
				data.type.int=FALSE, hessian=FALSE, unif.seed=optimcontrol$unif.seed, 
				int.seed=optimcontrol$int.seed,print.level=optimcontrol$print.level, share.type=0, instance.number=0,
				output.path="stdout", output.append=FALSE, project.path=NULL,
				P1=optimcontrol$P1, P2=optimcontrol$P2, P3=optimcontrol$P3, 
				P4=optimcontrol$P4, P5=optimcontrol$P5, P6=optimcontrol$P6,
				P7=optimcontrol$P7, P8=optimcontrol$P8, P9=optimcontrol$P9,
				P9mix=NULL, BFGSburnin=optimcontrol$BFGSburnin,BFGSfn=NULL, BFGShelp=NULL,
				cluster=FALSE, balance=FALSE, debug=FALSE,
				model=model, T=T, integration.points=integration.points,
				intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
				precalc.data=precalc.data,
				integration.weights=integration.weights,new.noise.var=new.noise.var,
				weight=weight)
				#hessian=TRUE
		
		o$par <- t(as.matrix(o$par))
		colnames(o$par) <- colnames(model@X)
		o$value <- as.matrix(o$value)
		colnames(o$value) <- colnames(model@y)
				
		return(list(par=o$par, value=o$value)) 
	}	
}
