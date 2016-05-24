summary.msgps <-
function(object, digits=max(3, getOption("digits") - 3), num.result = 20, coef.result=100,...){
		#check digits
		if(mode(digits)!="numeric")	stop('"digits" must be numeric.')
		if(length(digits) > 1)	stop('"digits" must be a scalar (1-dimensional vector).')
		if(as.integer(digits)!=digits)	stop('"digits" must be integer.')
		if(digits <= 0)	stop('"digits" must be positive integer.')

		#check num.result
		if(mode(num.result)!="numeric")	stop('"num.result" must be numeric.')
		if(length(num.result) > 1)	stop('"num.result" must be a scalar (1-dimensional vector).')
		if(as.integer(num.result)!=num.result)	stop('"num.result" must be integer.')
		if(num.result <= 0)	stop('"num.result" must be positive integer.')

		#check coef.result
		if(mode(coef.result)!="numeric")	stop('"coef.result" must be numeric.')
		if(length(coef.result) > 1)	stop('"coef.result" must be a scalar (1-dimensional vector).')
		if(as.integer(coef.result)!=coef.result)	stop('"coef.result" must be integer.')
		if(coef.result <= 0)	stop('"coef.result" must be positive integer.')
		
		
		STEP_adj <- object$dfgps_result$STEP_adj
		N_candidate <- min(num.result,STEP_adj)
		candidate_index <- trunc(seq(1,STEP_adj,length=N_candidate))
		result_tuning <- object$dfgps_result$tuning[candidate_index]
		if(object$stand.coef==FALSE) result_tuning <- object$dfgps_result$tuning_stand[candidate_index]
		result_df <- object$dfgps_result$df[candidate_index]

		result_df_bind <- cbind(result_tuning,result_df)
		colnames(result_df_bind) <- c("tuning","df")

		result_cp <-list(coef= object$dfcp_result$coef,tuning= object$dfcp_result$tuning,df= object$dfcp_result$df)
		result_aicc <-list(coef= object$dfaicc_result$coef,tuning= object$dfaicc_result$tuning,df= object$dfaicc_result$df)
		result_gcv <-list(coef= object$dfgcv_result$coef,tuning= object$dfgcv_result$tuning,df= object$dfgcv_result$df)
		result_bic <-list(coef= object$dfbic_result$coef,tuning= object$dfbic_result$tuning,df= object$dfbic_result$df)

		#coef.result <- object$coef.result
		p <- object$dfgps_result$p
		if(coef.result>=p) result_coef_ms <- cbind(object$dfcp_result$coef, object$dfaicc_result$coef, object$dfgcv_result$coef, object$dfbic_result$coef)
		result_tuning_ms <- cbind(object$dfcp_result$tuning, object$dfaicc_result$tuning, object$dfgcv_result$tuning, object$dfbic_result$tuning)
		result_df_ms <- cbind(object$dfcp_result$df, object$dfaicc_result$df, object$dfgcv_result$df, object$dfbic_result$df)
        if(coef.result>=p) 		colnames(result_coef_ms) <-  c("Cp","AICC","GCV","BIC")
		colnames(result_tuning_ms) <-  c("Cp","AICC","GCV","BIC")
		colnames(result_df_ms) <-  c("Cp","AICC","GCV","BIC")
		#tuning_maxの計算
		if(object$stand.coef==TRUE)  tuning.max = max(object$dfgps_result$tuning)
		if(object$stand.coef==FALSE)  tuning.max = max(object$dfgps_result$tuning_stand)
		
		penalty <- object$penalty
		alpha <- object$alpha
		gamma <- object$gamma
		lambda <- object$lambda

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
if(coef.result>=p) rownames(result_coef_ms) <- Rnames



#RESULT

   cat("\nCall:", paste(deparse(object$call)), "\n")
   cat("\nPenalty:", deparse(penalty), "\n")
    if(penalty == "enet" || penalty == "genet" ) cat("\nalpha:", format(alpha,digits=digits), "\n")
    if(penalty == "alasso" ) cat("\ngamma:", format(gamma,digits=digits), "\n")
    if(penalty == "alasso" ) cat("\nlambda:", format(lambda,digits=digits), "\n")
    cat("\ndf:\n"); print(result_df_bind,digits=digits);
    cat("\ntuning.max:",format(tuning.max,digits=digits),"\n")
    if(coef.result>=p) cat("\nms.coef:\n");
    if(coef.result>=p) print(result_coef_ms,digits=digits);
    cat("\nms.tuning:\n"); print(result_tuning_ms,digits=digits);
    cat("\nms.df:\n"); print(result_df_ms,digits=digits);
    cat("\n")
 
 	invisible(object)   
    
#		if(coef.result>=p){
#			if(penalty == "enet" || penalty == "genet" ) return(list(penalty=penalty, alpha=alpha, df=result_df_bind,tuning.max = tuning.max, ms.coef=result_coef_ms,ms.tuning=result_tuning_ms,ms.df=result_df_ms))
#			if(penalty == "alasso" ) return(list(penalty=penalty, gamma=gamma,lambda=lambda, df=result_df_bind,tuning.max = tuning.max, ms.coef=result_coef_ms,ms.tuning=result_tuning_ms,ms.df=result_df_ms))
#		}else{
#			if(penalty == "enet" || penalty == "genet" ) return(list(penalty=penalty, alpha=alpha, df=result_df_bind,tuning.max = tuning.max,ms.tuning=result_tuning_ms,ms.df=result_df_ms))
#			if(penalty == "alasso" ) return(list(penalty=penalty, gamma=gamma,lambda=lambda, df=result_df_bind,tuning.max = tuning.max,ms.tuning=result_tuning_ms,ms.df=result_df_ms))
#		}


}

