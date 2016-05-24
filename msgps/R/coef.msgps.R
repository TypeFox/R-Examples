coef.msgps <-
function(object,tuning,...){
	  #check tuning
	if(missing(tuning)==FALSE){
		if(mode(tuning)!="numeric")	stop('"tuning" must be numeric.')
		if(sum(tuning < 0)>=1)	stop('"tuning" must be non-negative.')
	}



	flag <- 0
	dfgps_result0 <- object$dfgps_result
	if(missing(tuning)){
		STEP <- rep(0,4)
		STEP[1] <- object$dfcp_result$step
		STEP[2] <- object$dfaicc_result$step
		STEP[3] <- object$dfgcv_result$step
		STEP[4] <- object$dfbic_result$step
		flag <- 1
	}
	if(flag==0) ans <- coef.dfgps(dfgps_result0, tuning, object$intercept, object$stand.coef)$coefficient
	if(flag==1) ans <- coef.step.dfgps(dfgps_result0, STEP, object$intercept, object$stand.coef)$coefficient
#	if(flag==1) names(ans$tuning) <- c("Cp","AICC","GCV","BIC")

#Rnames
	if(is.null(colnames(object$dfgps_result$X))==TRUE){
		if(object$intercept==TRUE){
			p0 <- object$dfgps_result$p+1
			Rnames <- rep(0,p0)
			Rnames[1] <- "(Intercept)"
			Rnames[2:p0] <- paste("V",1:(p0-1),sep="")
		}else{
			p00 <- object$dfgps_result$p
			Rnames <- rep(0,p00)
			Rnames[1:p00] <- paste("V",1:p00,sep="")
		}
	}else{
		if(object$intercept==TRUE){
			p0 <- object$dfgps_result$p+1
			Rnames <- rep(0,p0)
			Rnames[1] <- "(Intercept)"
			Rnames[2:p0] <- colnames(object$dfgps_result$X)
		}else{
			p00 <- object$dfgps_result$p
			Rnames <- colnames(object$dfgps_result$X)
		}
	}	
rownames(ans) <- Rnames

	
	
#Cnames
	if(flag==1){
		colnames(ans) <- c("Cp","AICc","GCV","BIC")
	}else{
		colnames(ans) <- paste("t=",tuning,sep="")
	}
	
	return(ans)
}

