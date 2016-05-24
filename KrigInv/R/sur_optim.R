sur_optim <- function(x, integration.points,integration.weights=NULL,
						intpoints.oldmean,intpoints.oldsd,precalc.data,
						model, T, new.noise.var=NULL,current.sur=1e6){
	
	#integration.points = NULL : remplir avec integration_param				
					
	tp1<-as.numeric(t(model@X))
	tp2<-matrix(tp1-as.numeric(x),nrow=model@n,byrow=TRUE)^2
	mindist<-min(sqrt(rowSums(tp2)))
	if (!identical(colnames(integration.points), colnames(model@X))) colnames(integration.points) <- colnames(model@X)
	
	if ( (mindist > 1e-5) || ( !is.null(new.noise.var)&&(new.noise.var!=0)  ) ){
		X.new <- t(x)
		krig  <- predict_nobias_km(object=model, newdata=as.data.frame(X.new), 
									type="UK",se.compute=TRUE, 
									cov.compute=FALSE) 
		
		mk <- krig$mean ; sk <- krig$sd ; newXvar <- sk*sk
		F.newdata <- krig$F.newdata ; c.newdata <- krig$c
		kn = computeQuickKrigcov(model,integration.points,X.new,precalc.data, F.newdata , c.newdata) 
		
		krig  <- predict_update_km (newXmean=mk,newXvar=newXvar,newXvalue=mk, 
				newdata.oldmean=intpoints.oldmean,newdata.oldsd=intpoints.oldsd,kn=kn)
			
		sk.new <- krig$sd;lambda <- krig$lambda	
		a <- (intpoints.oldmean-T) / sk.new
		b <- lambda * sk / sk.new
		
		a[a==Inf]<- 1000 ;a[a== -Inf] <- 1000;b[b==Inf]<-0;b[b== -Inf]<-0
		a[is.nan(a)] <- 1000; b[is.nan(b)]<-0
		bsquare <- b*b
		a.new <- as.numeric(a/sqrt(1+bsquare))
		b.new <- as.numeric(-bsquare/(1+bsquare))
		Phi.biv.a.b <- pbivnorm(a.new,-a.new,b.new)  #c.d.f of the bivariate gaussian distribution
		
		if (is.null(integration.weights)) {crit <- mean(Phi.biv.a.b)
		}else crit <- sum(Phi.biv.a.b*integration.weights)
	}else crit <- current.sur		
	
	return(crit)
}
