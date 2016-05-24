timse_optim <- function(x, integration.points,integration.weights=NULL,
		intpoints.oldmean=NULL,intpoints.oldsd=NULL,precalc.data,
		model, T, new.noise.var=0,weight=NULL){
	
	if(!is.null(new.noise.var)){
		if(new.noise.var==0) {
			new.noise.var <- NULL
		}
	}					
	
	tp1<-as.numeric(t(model@X))
	tp2<-matrix(tp1-as.numeric(x),nrow=model@n,byrow=TRUE)^2
	mindist<-min(sqrt(rowSums(tp2)))
	if (!identical(colnames(integration.points), colnames(model@X))) colnames(integration.points) <- colnames(model@X)
	
	if ((mindist > 1e-5) || (!is.null(new.noise.var))){
		X.new <- t(x)
		krig  <- predict_nobias_km(object=model, newdata=as.data.frame(X.new), 
								type="UK",se.compute=TRUE,cov.compute=FALSE) 
		
		mk <- krig$mean;sk <- krig$sd;newXvar<-sk*sk
		F.newdata<-krig$F.newdata ; c.newdata<-krig$c
		kn = computeQuickKrigcov(model,integration.points,X.new,precalc.data, F.newdata , c.newdata) 
		
		krig  <- predict_update_km (newXmean=mk,newXvar=newXvar,newXvalue=mk, 
				newdata.oldmean=intpoints.oldmean,newdata.oldsd=intpoints.oldsd,kn=kn)
		
		sk.new <- krig$sd	
		if(is.null(weight)){
			timse <- sk.new^2
		}else{
			timse <- weight * sk.new^2
		}
		
		if (is.null(integration.weights)) {crit <- mean(timse)
		}else crit <- sum(timse*integration.weights)
	
	}else crit <- Inf		
	
	return(crit)
	
}
	
