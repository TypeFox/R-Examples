timse_optim_parallel <- function(x, integration.points,integration.weights=NULL,
		intpoints.oldmean=NULL,intpoints.oldsd=NULL,precalc.data,
		model, T, new.noise.var=0,weight=NULL,batchsize,current.timse){
	
	if(!is.null(new.noise.var)){
		if(new.noise.var==0) {
			new.noise.var <- NULL
		}
	}					
	
	#x is a vector of size d * batchsize
	d <- model@d
	n <- model@n
	X.new <- matrix(x,nrow=d)
	mindist <- Inf

	tp1 <- c(as.numeric(t(model@X)),x)
	for (i in 1:batchsize){
		#distance between the i^th point and all other points (in the DOE or in the batch)
		xx <- X.new[,i]
		tp2<-matrix(tp1-as.numeric(xx),ncol=d,byrow=TRUE)^2
		mysums <- sqrt(rowSums(tp2))
		mysums[n+i] <- Inf #because this one is usually equal to zero...
		mindist <- min(mindist,mysums)		
	}
		
	if (!identical(colnames(integration.points), colnames(model@X))) colnames(integration.points) <- colnames(model@X)
	
	if ((mindist > 1e-5) || (!is.null(new.noise.var))){
		X.new <- t(X.new)	
		krig  <- predict_nobias_km(object=model, newdata=as.data.frame(X.new), 
								type="UK",se.compute=TRUE, cov.compute=TRUE) 
	
		mk <- krig$mean ; sk <- krig$sd ; newXvar <- sk*sk
		F.newdata <- krig$F.newdata ; c.newdata <- krig$c;Sigma.r <- krig$cov
		Sigma.r.norm <- krig$cov / crossprod(t(sk))
	
		kn = computeQuickKrigcov(model,integration.points,X.new,precalc.data, F.newdata , c.newdata) 
	
		krig2  <- predict_update_km_parallel (newXmean=mk,newXvar=newXvar,newXvalue=mk, 
				Sigma.r=Sigma.r,newdata.oldmean=intpoints.oldmean,newdata.oldsd=intpoints.oldsd,kn=kn)
		if(!is.null(krig2$error)) return(current.timse)

		sk.new <- krig2$sd	
		if(is.null(weight)){
			tmse <- sk.new^2	
		}else{
			tmse <- weight * sk.new^2
		}

		if (is.null(integration.weights)) {crit <- mean(tmse)
		}else crit <- sum(tmse*integration.weights)

	}else crit <- current.timse + 0.01		

	return(crit)
	
}
