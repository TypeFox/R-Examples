 max_vorob_parallel <- function(lower, upper, optimcontrol=NULL, 
		 				batchsize,
						integration.param,
						T, model, new.noise.var=0
            ){
	
  optim.option <- optimcontrol$optim.option
  if(is.null(optim.option)) optim.option <- 2
  
	integration.points <- as.matrix(integration.param$integration.points) ; d <- model@d
	integration.weights <- integration.param$integration.weights
  alpha <- integration.param$alpha
	
  if(is.null(optimcontrol$method)) optimcontrol$method <- "genoud"
	
	#precalculates the kriging mean and variance on the integration points
	pred <- predict_nobias_km(object=model,newdata=integration.points,type="UK",se.compute=TRUE)
	intpoints.oldmean <- pred$mean ; intpoints.oldsd <- pred$sd; pn <- pnorm((intpoints.oldmean-T)/intpoints.oldsd)
  
  if(is.null(alpha)) alpha <- vorob_threshold(pn)
  pn_bigger_than_alpha <- (pn>alpha)+0
  pn_lower_than_alpha <- 1-pn_bigger_than_alpha
  
	if(is.null(integration.weights)) current.vorob <- mean(pn*pn_lower_than_alpha + (1-pn)*pn_bigger_than_alpha)
	if(!is.null(integration.weights)) current.vorob <- sum(integration.weights*(pn*pn_lower_than_alpha + (1-pn)*pn_bigger_than_alpha))
	precalc.data <- precomputeUpdateData(model,integration.points)
				
	fun.optim <- vorob_optim_parallel
	########################################################################################
	#discrete Optimisation
	#batchsize optimizations in dimension d
	if(optimcontrol$method=="discrete"){
		
		if (is.null(optimcontrol$optim.points)){
			n.discrete.points <- d*100
			optimcontrol$optim.points <- t(lower + t(matrix(runif(d*n.discrete.points),ncol=d)) * (upper - lower))
		}
    
		optim.points <- optimcontrol$optim.points
		optim.points <- data.frame(optim.points)
    
		if(ncol(optim.points)==d){
		  #this is the standard case:
		  fun.optim <- vorob_optim_parallel2
		  colnames(optim.points) <- colnames(model@X)
		  all.crit <- seq(1,nrow(optim.points))
		  if(nrow(optim.points) < batchsize){
			  print("error in max_vorob_parallel")
			  print("please set a batchsize lower or equal than the number of tested points optimcontrol$optim.points")
		  }
		
		  other.points <- NULL
		  for (j in 1:batchsize){
						
			  for (i in 1:nrow(optim.points)){
				  all.crit[i] <- fun.optim(x=t(optim.points[i,]), integration.points=integration.points,integration.weights=integration.weights,
					intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
					precalc.data=precalc.data,T=T, model=model, new.noise.var=new.noise.var,
					other.points=other.points,batchsize=j,alpha=alpha,current.vorob=current.vorob)
			  }	
			  ibest <- which.min(all.crit)
			  other.points <- c(other.points,as.numeric(optim.points[ibest,]))					
		  }
			
		  o <- list(3)
		  o$par <- other.points;o$par <- t(matrix(o$par,nrow=d)); colnames(o$par) <- colnames(model@X)
		  o$value <- min(all.crit); o$value <- as.matrix(o$value); colnames(o$value) <- colnames(model@y)
		  o$allvalues <- all.crit
		  return(list(par=o$par, value=o$value,allvalues=o$allvalues))
		}else{
		  #new code (Aug 2012)
		  fun.optim <- vorob_optim_parallel2
		  all.crit <- seq(1,nrow(optim.points))
      
		  for (i in 1:nrow(optim.points)){
		    all.crit[i] <- fun.optim(x=t(optim.points[i,]), integration.points=integration.points,integration.weights=integration.weights,
		                             intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
		                             precalc.data=precalc.data,T=T, model=model, new.noise.var=new.noise.var,
		                             batchsize=batchsize,alpha=alpha,current.vorob=current.vorob)
		  }
		  ibest <- which.min(all.crit)
		  o <- list(3)
		  o$par <- t(matrix(optim.points[ibest,],nrow=d,ncol=batchsize)); colnames(o$par) <- colnames(model@X)
		  o$value <- all.crit[ibest]; o$value <- as.matrix(o$value); colnames(o$value) <- colnames(model@y)
		  o$allvalues <- all.crit
		  return(list(par=o$par, value=o$value,allvalues=o$allvalues))
		}
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
		if (is.null(optimcontrol$P1)) optimcontrol$P1<-0#50
		if (is.null(optimcontrol$P2)) optimcontrol$P2<-0#50
		if (is.null(optimcontrol$P3)) optimcontrol$P3<-0#50
		if (is.null(optimcontrol$P4)) optimcontrol$P4<-0#50
		if (is.null(optimcontrol$P5)) optimcontrol$P5<-50
		if (is.null(optimcontrol$P6)) optimcontrol$P6<-50#50
		if (is.null(optimcontrol$P7)) optimcontrol$P7<-50
		if (is.null(optimcontrol$P8)) optimcontrol$P8<-50
		if (is.null(optimcontrol$P9)) optimcontrol$P9<-0
		
		if(optim.option==1){
			#one unique optimisation in dimension batchsize * d 
			domaine <- cbind(rep(lower,times=batchsize), rep(upper,times=batchsize))
			
			o <- genoud(fn=fun.optim, nvars=d*batchsize, max=FALSE, pop.size=optimcontrol$pop.size,
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
				precalc.data=precalc.data,integration.weights=integration.weights,
				new.noise.var=new.noise.var,batchsize=batchsize,alpha=alpha,
        current.vorob=current.vorob)
			
			o$par <- t(matrix(o$par,nrow=d)); colnames(o$par) <- colnames(model@X)
			o$value <- as.matrix(o$value); colnames(o$value) <- colnames(model@y)
		}else{
			#batchsize optimisations in dimension d
			fun.optim <- vorob_optim_parallel2
			domaine <- cbind(lower,upper)
			other.points <- NULL
			for (i in 1:batchsize){
		
				o <- genoud(fn=fun.optim, nvars=d, max=FALSE, pop.size=optimcontrol$pop.size,
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
						cluster=FALSE, balance=FALSE, debug=FALSE,other.points=other.points,
						model=model, T=T, integration.points=integration.points,
						intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
						precalc.data=precalc.data,integration.weights=integration.weights,
						new.noise.var=new.noise.var,batchsize=i,current.vorob=current.vorob,
            alpha=alpha)
        
        other.points <- c(other.points,as.numeric(o$par))
			}
			o$par <- t(matrix(other.points,nrow=d)); colnames(o$par) <- colnames(model@X)
			o$value <- as.matrix(o$value); colnames(o$value) <- colnames(model@y)
		}
				
		return(list(par=o$par, value=o$value)) 
	}	
}

