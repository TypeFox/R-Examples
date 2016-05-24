plot.df <-function(object, criterion="cp",  stand.coef=TRUE, plot.step = 1000, col=TRUE){		if(criterion!="cp" && criterion!="aicc" && criterion!="gcv" && criterion!="bic" && criterion!="none")	stop('"criterion" must be "cp", "aicc", "gcv", "bic" or "none".')		if(!is.numeric(plot.step))	stop("plot.step must be a numeric.")		if(length(plot.step)>1)	stop("plot.step must be a scalar.")		if(plot.step<5)	stop("plot.step must be greater than 5.")		if(col!=TRUE && col!=FALSE) stop('"col" must be TRUE or FALSE')		if(stand.coef!=TRUE && stand.coef!=FALSE) stop('"col" must be TRUE or FALSE')#        if(sum(tuning<0)>0)	stop('"tuning" must be non-negative.')#		if(stand.coef!=TRUE && stand.coef!=FALSE) stop('"stand.coef" must be TRUE or FALSE')#		if(intercept!=TRUE && intercept!=FALSE) stop('"intercept" must be TRUE or FALSE')		STEP_adj <- object$dfgps_result$STEP_adj		N_candidate <- min(plot.step,STEP_adj)		candidate_index <- trunc(seq(1,STEP_adj,length=N_candidate))
		dfgps_result0 <- object$dfgps_result
		coefficient <- coef.step.dfgps(dfgps_result0, candidate_index, intercept=FALSE, stand.coef= stand.coef)
		tuning2 <- apply(abs(coefficient$coefficient),2,sum)
		maxtuning2 <- max(tuning2)
		tuning2 <- tuning2 / maxtuning2
		#result_tuning <- object$dfgps_result$tuning[candidate_index]		#if(stand.coef==FALSE) result_tuning <- object$dfgps_result$tuning_stand[candidate_index]		result_df <- object$dfgps_result$df[candidate_index]		#result_df_bind <- cbind(result_tuning,result_df)
		result_df_bind <- cbind(tuning2,result_df)
		xrange <- c(0,1)
		yrange <- c(0,object$dfgps_result$p)
		plot(result_df_bind,xlab="|beta|/max|beta|",ylab="df",type="l",xlim= xrange,ylim= yrange)
		if(object$alpha==0 &&  object$penalty=="enet"){
			dfZou <- apply(abs(coefficient$coefficient) > object$dfgps_result$delta_t*1.5 ,2,sum)
			par(new=T)
			plot(tuning2,dfZou,ann=F,xlab=",ylab=",type="l",xlim= xrange,ylim= yrange,lty=2)
		}
		
		
		}