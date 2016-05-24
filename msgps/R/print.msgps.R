print.msgps=function(x,digits = max(3, getOption("digits") - 3),num.result = 20,...){
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

		STEP_adj <- x$dfgps_result$STEP_adj
		N_candidate <- min(num.result,STEP_adj)
		candidate_index <- trunc(seq(1,STEP_adj,length=N_candidate))
		result_tuning <- x$dfgps_result$tuning[candidate_index]
		if(x$stand.coef==FALSE) result_tuning <- x$dfgps_result$tuning_stand[candidate_index]
		result_df <- x$dfgps_result$df[candidate_index]

		result_df_bind <- cbind(result_tuning,result_df)
		colnames(result_df_bind) <- c("tuning","df")

		result_cp <-list(coef= x$dfcp_result$coef,tuning= x$dfcp_result$tuning,df= x$dfcp_result$df)
		result_aicc <-list(coef= x$dfaicc_result$coef,tuning= x$dfaicc_result$tuning,df= x$dfaicc_result$df)
		result_gcv <-list(coef= x$dfgcv_result$coef,tuning= x$dfgcv_result$tuning,df= x$dfgcv_result$df)
		result_bic <-list(coef= x$dfbic_result$coef,tuning= x$dfbic_result$tuning,df= x$dfbic_result$df)

		#coef.result <- x$coef.result
		#p <- x$dfgps_result$p
		#if(coef.result>=p) result_coef_ms <- cbind(x$dfcp_result$coef, x$dfaicc_result$coef, x$dfgcv_result$coef, x$dfbic_result$coef)
		result_tuning_ms <- cbind(x$dfcp_result$tuning, x$dfaicc_result$tuning, x$dfgcv_result$tuning, x$dfbic_result$tuning)
		result_df_ms <- cbind(x$dfcp_result$df, x$dfaicc_result$df, x$dfgcv_result$df, x$dfbic_result$df)
		colnames(result_tuning_ms) <-  c("Cp","AICC","GCV","BIC")
		colnames(result_df_ms) <-  c("Cp","AICC","GCV","BIC")
		#tuning_maxの計算
		if(x$stand.coef==TRUE)  tuning.max = max(x$dfgps_result$tuning)
		if(x$stand.coef==FALSE)  tuning.max = max(x$dfgps_result$tuning_stand)
		
		penalty <- x$penalty
		alpha <- x$alpha
		gamma <- x$gamma
		lambda <- x$lambda

#RESULT

   cat("\nCall:", paste(deparse(x$call)), "\n")
   cat("\nPenalty:", deparse(penalty), "\n")
    if(penalty == "enet" || penalty == "genet" ) cat("\nalpha:", format(alpha,digits=digits), "\n")
    if(penalty == "alasso" ) cat("\ngamma:", format(gamma,digits=digits), "\n")
    if(penalty == "alasso" ) cat("\nlambda:", format(lambda,digits=digits), "\n")
    #cat("\ndf:\n"); print(result_df_bind,digits=digits);
    cat("\nms.tuning:\n"); print(result_tuning_ms,digits=digits);
    cat("\n")
 
 	invisible(x)   
    
#		if(coef.result>=p){
#			if(penalty == "enet" || penalty == "genet" ) return(list(penalty=penalty, alpha=alpha, df=result_df_bind,tuning.max = tuning.max, ms.coef=result_coef_ms,ms.tuning=result_tuning_ms,ms.df=result_df_ms))
#			if(penalty == "alasso" ) return(list(penalty=penalty, gamma=gamma,lambda=lambda, df=result_df_bind,tuning.max = tuning.max, ms.coef=result_coef_ms,ms.tuning=result_tuning_ms,ms.df=result_df_ms))
#		}else{
#			if(penalty == "enet" || penalty == "genet" ) return(list(penalty=penalty, alpha=alpha, df=result_df_bind,tuning.max = tuning.max,ms.tuning=result_tuning_ms,ms.df=result_df_ms))
#			if(penalty == "alasso" ) return(list(penalty=penalty, gamma=gamma,lambda=lambda, df=result_df_bind,tuning.max = tuning.max,ms.tuning=result_tuning_ms,ms.df=result_df_ms))
#		}


}
