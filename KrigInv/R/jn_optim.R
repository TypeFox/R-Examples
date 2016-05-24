jn_optim <- function(x, integration.points,integration.weights=NULL,
					intpoints.oldmean,intpoints.oldsd,precalc.data,
					model, T, new.noise.var=NULL,current.sur=0){
	
					
	
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
		if(length(integration.points) == (model@d * 2 * length(integration.weights))){
			#special case where integration points are chosen with a jn distribution
			M <- nrow(integration.points) / 2
			
			a <- (intpoints.oldmean-T) / sk.new
			b <- lambda * sk / sk.new
			a[a==Inf]<- 1000 ;a[a== -Inf] <- 1000;b[b==Inf]<-0;b[b== -Inf]<-0
			a[is.nan(a)] <- 1000; b[is.nan(b)]<-0
			
			a1 <- a[1:M] ; a2 <- a[(M+1):(2*M)]
			b1 <- b[1:M] ; b2 <- b[(M+1):(2*M)]
			tmp1 <- sqrt(1+b1*b1) ; tmp2 <- sqrt(1+b2*b2)
			arg1 <-  as.numeric(a1/tmp1)
			arg2 <-  as.numeric(a2/tmp2)
			arg3 <- as.numeric(b1*b2/(tmp1*tmp2))
			arg3 <- pmin(arg3,1);arg3 <- pmax(arg3,-1)
			Tau.a.b <- pbivnorm(arg1,arg2,arg3)
			
			crit <- sum(Tau.a.b*integration.weights) * (-1)
			
		}else{
			#general case with M^2 integration points 
			a <- (intpoints.oldmean-T) / sk.new
			b <- lambda * sk / sk.new
			
			a[a==Inf]<- 1000 ;a[a== -Inf] <- 1000;b[b==Inf]<-0;b[b== -Inf]<-0
			a[is.nan(a)] <- 1000; b[is.nan(b)]<-0
			bsquare <- b*b
			M <- nrow(integration.points)
			indices <- expand.grid(c(1:M),c(1:M))
		
			col1 <- indices[,1] ; col2 <- indices[,2]
			tmp1 <- sqrt(1+bsquare[col1]) ; tmp2 <- sqrt(1+bsquare[col2])
		
			arg1 <-  as.numeric(a[col1]/tmp1)
			arg2 <-  as.numeric(a[col2]/tmp2)
			arg3 <- as.numeric(b[col1]*b[col2]/(tmp1*tmp2))
      arg3 <- pmin(arg3,1);arg3 <- pmax(arg3,-1)
			Tau.a.b <- pbivnorm(arg1,arg2,arg3)
					
			if (is.null(integration.weights)) {crit <- mean(Tau.a.b)* (-1)
			}else{ 
				new.weights <- integration.weights[col1]*integration.weights[col2]
				crit <- sum(Tau.a.b*new.weights) * (-1)
			}
		}
	}else crit <- current.sur		
	
	return(crit)
}
