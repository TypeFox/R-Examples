plot.msgps <-
function(x, criterion="cp",  xvar="norm", yvar="coef", yvar.dflasso=TRUE, stand.coef=TRUE, plot.step = 1000, col=TRUE,...){
		if(criterion!="cp" && criterion!="aicc" && criterion!="gcv" && criterion!="bic" && criterion!="none")	stop('"criterion" must be "cp", "aicc", "gcv", "bic" or "none".')
		if(xvar!="norm" && xvar!="sum" && xvar!="step" &&  xvar!="t") stop('"xvar" must be "norm", "sum", "step", or "t".')
		if(yvar!="coef" && yvar!="df") stop('"yvar" must be "coef", "df".')
		if(!is.numeric(plot.step))	stop("plot.step must be a numeric.")
		if(length(plot.step)>1)	stop("plot.step must be a scalar.")
		if(plot.step<5)	stop("plot.step must be greater than 5.")
		if(col!=TRUE && col!=FALSE) stop('"col" must be TRUE or FALSE')
		if(stand.coef!=TRUE && stand.coef!=FALSE) stop('"stand.coef" must be TRUE or FALSE')
		if(yvar.dflasso!=TRUE && yvar.dflasso!=FALSE) stop('"yvar.dflasso" must be TRUE or FALSE')
		if(mode(plot.step)!="numeric")	stop('"plot.step" must be numeric.')
		if(length(plot.step) > 1)	stop('"plot.step" must be a scalar (1-dimensional vector).')
		if(as.integer(plot.step)!=plot.step)	stop('"plot.step" must be integer.')
		if(plot.step <= 0)	stop('"plot.step" must be positive integer.')
		if(plot.step >= 1e+6)	stop('"plot.step" must be less than 1e+6.')

#        if(sum(tuning<0)>0)	stop('"tuning" must be non-negative.')
#		if(stand.coef!=TRUE && stand.coef!=FALSE) stop('"stand.coef" must be TRUE or FALSE')
#		if(intercept!=TRUE && intercept!=FALSE) stop('"intercept" must be TRUE or FALSE')

			#tuning.max <- max(x$dfgps_result$tuning)
			#if(stand.coef == FALSE) tuning.max  <- max(x$dfgps_result$tuning_stand)
			#tuning <- seq(0,tuning.max,length=plot.step)
			#dfgps_result0 <- x$dfgps_result
			#coefficient <- coef.dfgps(dfgps_result0, tuning, intercept=FALSE, stand.coef= stand.coef)
			#coefficient <- coef(x,tuning)
			#tuning2 <- coefficient$tuning
			#tuning2_rank <- rank(tuning2,ties.method="first")

	
			STEPmax <- x$dfgps_result$STEP_adj
	#		if(stand.coef == FALSE) tuning.max  <- max(x$dfgps_result$tuning_stand)
			N_candidate <- min(plot.step, STEPmax)
			STEP_vec <- trunc(seq(1,STEPmax,length=N_candidate))
			tuning <- x$dfgps_result$tuning[STEP_vec]
			dfgps_result0 <- x$dfgps_result
			coefficient <- coef.step.dfgps(dfgps_result0, STEP_vec, intercept=FALSE, stand.coef= stand.coef)
	#		coefficient <- coef(x,tuning)
	#		tuning2 <- coefficient$tuning
			tuning2 <- apply(abs(coefficient$coefficient),2,sum)
	
			if(xvar=="norm"){
				maxtuning2 <- max(tuning2)
				tuning2 <- tuning2 / maxtuning2
				}
			
			if(xvar=="step") tuning2 <- STEP_vec
			if(xvar=="t") tuning2 <- tuning

#		tuning2_rank <- rank(tuning2,ties.method="first")


	if(xvar=="norm") xlab="|beta|/max|beta|"
	if(xvar=="sum") xlab="|beta|"
	if(xvar=="step") xlab="step"
	if(xvar=="t") xlab="t"


####STEP_vec修正

	if(xvar!="step"){
		tuning2_rank <- rank(tuning2,ties.method="first")
		tuning2 <- tuning2[tuning2_rank]
		tuning2 <- tuning2[tuning2_rank]
		STEP_vec <- STEP_vec[tuning2_rank]
	}


#############

########df
	if(yvar=="df"){
		result_df <- x$dfgps_result$df[STEP_vec]
		#result_df_bind <- cbind(result_tuning,result_df)
		result_df_bind <- cbind(tuning2,result_df)
		yrange <- c(0,x$dfgps_result$p)
		plot(result_df_bind,xlab=xlab, ylab="df",type="l", ylim= yrange,...)
		if(x$alpha==0 &&  x$penalty=="enet" && yvar.dflasso==TRUE){
			dfZou <- apply(abs(coefficient$coefficient) > x$dfgps_result$delta_t*1.5 ,2,sum)
			par(new=T)
			plot(tuning2,dfZou,ann=F,xlab="",ylab="",type="l",ylim= yrange,axes = F, lty=2,...)
		}
	}
###########



###########solution path
	if(yvar=="coef"){
		#		tuning2 <- apply(abs(coefficient$coefficient),2,sum)
				if(stand.coef==TRUE){
					if(col==TRUE) plot(tuning2,coefficient$coefficient[1,],ylim=range(coefficient$coefficient),xlab=xlab, ylab="standardized coefficients",type="l",lty=1, col=1,...)
					if(col==FALSE) plot(tuning2,coefficient$coefficient[1,],ylim=range(coefficient$coefficient), xlab=xlab, ylab="standardized coefficients", type="l",  lty=1,...)
				}
				if(stand.coef==FALSE){
					if(col==TRUE) plot(tuning2,coefficient$coefficient[1,],ylim=range(coefficient$coefficient),xlab=xlab, ylab="coefficients",type="l",lty=1, col=1,...)
					if(col==FALSE) plot(tuning2,coefficient$coefficient[1,],ylim=range(coefficient$coefficient), xlab=xlab, ylab="coefficients", type="l",  lty=1,...)
				}
				par(new=T)
				for(i in 2:nrow(coefficient$coefficient)){
					if(col==TRUE) plot(tuning2,coefficient$coefficient[i,],ylim=range(coefficient$coefficient),ann=F,type="l", axes = F,lty=i+1, col=i+1,...)
					if(col==FALSE) plot(tuning2,coefficient$coefficient[i,],ylim=range(coefficient$coefficient),ann=F,type="l", axes = F, lty=i+1,...)
					par(new=T)
			}
			par(new=F)
			abline(h=0,lwd=1,col="black")
			if(criterion=="cp"){
				v00 <- coef.step.dfgps(dfgps_result0, x$dfcp_result$step, intercept=FALSE, stand.coef= stand.coef)
				if(xvar=="step") v0 <- x$dfcp_result$step
				#if(col==TRUE) abline(v=x$dfgps_result$tuning[x$dfcp_result$step],col="red")
				#if(col==FALSE) abline(v=x$dfgps_result$tuning[x$dfcp_result$step])
			}
			if(criterion=="aicc"){
				v00 <- coef.step.dfgps(dfgps_result0, x$dfaicc_result$step, intercept=FALSE, stand.coef= stand.coef)
				if(xvar=="step") v0 <- x$dfaicc_result$step
				#if(col==TRUE) abline(v=x$dfgps_result$tuning[x$dfaicc_result$step],col="red")
				#if(col==FALSE) abline(v=x$dfgps_result$tuning[x$dfaicc_result$step])
			}
			if(criterion=="gcv"){
				v00 <- coef.step.dfgps(dfgps_result0, x$dfgcv_result$step, intercept=FALSE, stand.coef= stand.coef)
				if(xvar=="step") v0 <- x$dfgcv_result$step
				#if(col==TRUE) abline(v=x$dfgps_result$tuning[x$dfgcv_result$step],col="red")
				#if(col==FALSE) abline(v=x$dfgps_result$tuning[x$dfgcv_result$step])
			}
			if(criterion=="bic"){
				v00 <- coef.step.dfgps(dfgps_result0, x$dfbic_result$step, intercept=FALSE, stand.coef= stand.coef)
				if(xvar=="step") v0 <- x$dfbic_result$step
				#if(col==TRUE) abline(v=x$dfgps_result$tuning[x$dfbic_result$step],col="red")
				#if(col==FALSE) abline(v=x$dfgps_result$tuning[x$dfbic_result$step])
			}
		
		
			if(xvar=="norm" && criterion!="none") v0 <- sum(abs(v00$coefficient)) / maxtuning2
			if(xvar=="sum" && criterion!="none") v0 <- sum(abs(v00$coefficient))
			if(xvar=="t" && criterion!="none") v0 <- v00$tuning
		
		
			if(criterion!="none"){
				if(col==TRUE) abline(v=v0,col="red",...)
				if(col==FALSE) abline(v=v0,...)
			}
	}
		
	#	if(type=="df"){ 
	#	STEP_adj <- x$dfgps_result$STEP_adj
	#	N_candidate <- min(plot.step,STEP_adj)
	#	candidate_index <- trunc(seq(1,STEP_adj,length=N_candidate))
	#	result_tuning <- x$dfgps_result$tuning[candidate_index]
	#	if(stand.coef==FALSE) result_tuning <- x$dfgps_result$tuning_stand[candidate_index]
	#	result_df <- x$dfgps_result$df[candidate_index]
	#	result_df_bind <- cbind(result_tuning,result_df)
	#	plot(result_df_bind,xlab="t",ylab="df",type="l")
	#	}
	
}