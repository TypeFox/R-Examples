coef.step.dfgps <-
function(object, STEP, intercept=TRUE, stand.coef=FALSE,...){
        #library.dynam("msgps")
		if(mode(STEP)!="numeric")	stop("STEP must be numeric.")
		if(!is.vector(STEP))	stop("STEP must be a vector.")
        if(sum(STEP<=0)>0)	stop('"STEP" must be a positive integer.')
		if(sum((STEP-1)>object$STEP_adj) >0) stop('"STEP" must not exceed the number of maximum step number.') #changed 2011.11.1
		if(stand.coef!=TRUE && stand.coef!=FALSE) stop('"stand.coef" must be TRUE or FALSE')
		if(intercept!=TRUE && intercept!=FALSE) stop('"intercept" must be TRUE or FALSE')

		tuning_candidate <- object$tuning
		if(stand.coef==FALSE) tuning_candidate<- object$tuning_stand
#		n_tuning <- length(tuning)
 		#tuning <- sort(tuning)
#		tuning <- as.numeric(tuning)
	    step0=STEP
#	   tuning_candidate_mat <- matrix(rep(tuning_candidate, n_tuning), n_tuning,length(tuning_candidate),byrow=T)
		#step0 <- apply(abs(tuning_candidate_mat-tuning),1, which.min)
#	   step0 <- sort(step0)
		step00 <- step0[step0==1]
		step01 <- step0[step0!=1]
		ans <- matrix(0,object$p,length(step0))

       if(length(step00)!=0){
                ans[,1:length(step00)]=rep(0,object$p)
       }
	   if(length(step01)!=0){
                step01 = as.integer(step01-1)
				step01_adj <- rep(0,max(step01))
				for(i in 1:length(step01)) step01_adj[step01[i]] <- step01_adj[step01[i]]+1
				step01_adj <- as.integer(step01_adj)
#			step01_order <- order(step01)
				step01_rank <- rank(step01,ties.method="first")
#       			step01_order <-  as.integer(step01_order-1)
                step01 = sort(step01)
                step01 = as.integer(step01)
				stepmax <- max(step01)
				stepmax <- as.integer(stepmax)
                 p = object$p
                 p = as.integer(p)
                delta_t=object$delta_t
                delta_t=as.double(delta_t)
                betamat = as.integer(object$coefficient_index)
                 beta0=.Call("betaOUT",betamat,step01,p,delta_t,step01_adj,stepmax)
				step01 <-step01+ 1
				beta0 <- beta0[,step01_rank]
				ans[,(length(step00)+1):length(step0)] <- beta0
				 }
		
			if(stand.coef==FALSE){
				temp <- sqrt(object$N-1) * sqrt(apply(object$X,2,var))
				temp_mat <- matrix(rep(temp,length(step0)),object$p,length(step0))
				ans <- ans / temp_mat
				ans[is.nan(ans)] <- 0
			}
			beta00 <- - t(apply(object$X,2,mean)) %*% ans + mean(object$y)
			if(intercept==TRUE) ans <- rbind(beta00,ans)
			

	
	
	
			 return(list(coefficient=ans,tuning=tuning_candidate[step0]))
}

