predict.msgps <- function(object,X,tuning,...){
	  #check X
	   if(!is.matrix(X) && !is.vector(X)) stop(" X must be a matrix of vector.")
   	 if(mode(X)!="numeric") stop(" X must be numeric.")
	   if (sum(complete.cases(X)==FALSE)>0)  stop("X must be complete data.")
	   if (is.matrix(X)==TRUE) if(ncol(X) != object$dfgps_result$p) stop("The number of variables is not correct.")
	   if(is.vector(X)==TRUE){
	   		if(length(X) != object$dfgps_result$p) stop("The number of variables is not correct.")
			X <- as.matrix(X)
			X <- t(X)
	   	}
		   
	  #check tuning
	if(missing(tuning)==FALSE){
		if(mode(tuning)!="numeric")	stop('"tuning" must be numeric.')
		if(sum(tuning < 0)>=1)	stop('"tuning" must be non-negative.')
	}	   	


if(missing(tuning)){
	if(object$intercept==TRUE){
			aicc <- object$dfaicc_result$coef[1] + X %*% object$dfaicc_result$coef[-1]
			bic <- object$dfbic_result $coef[1] + X %*% object$dfbic_result $coef[-1]
			gcv <- object$dfgcv_result $coef[1] + X %*% object$dfgcv_result $coef[-1]
			cp <- object$dfcp_result $coef[1] + X %*% object$dfcp_result $coef[-1]
	}
		if(object$intercept==FALSE){
			aicc <- X %*% object$dfaicc_result$coef
			bic <-  X %*% object$dfbic_result $coef
			gcv <- X %*% object$dfgcv_result $coef
			cp <- X %*% object$dfcp_result $coef
		}
		ans <- cbind(cp,aicc,gcv,bic)
		colnames(ans) <-  c("Cp","AICc","GCV","BIC")
	}else{
		 if(mode(tuning)!="numeric") stop(" tuning must be numeric.")
		 if(!is.vector(tuning)) stop(" tuning must be a vector.")
	 	if(sum(tuning<0)<0) stop(" tuning must be non-negative.")
			coef <-  coef.msgps(object,tuning)
			if(object$intercept==TRUE)  ans <- matrix(rep(coef[1,],nrow(X)),nrow(X),length(tuning),byrow=T) + X %*% coef[-1,]
			if(object$intercept==FALSE)  ans <- X %*% coef
	
		#Cnames
		colnames(ans) <- paste("t=",tuning,sep="")
		
	}
	
	
	return(ans)
}