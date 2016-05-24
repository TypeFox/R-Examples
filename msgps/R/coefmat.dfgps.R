coefmat.dfgps <-
function(dfgps_result, MAX.p = 200, intercept=TRUE, stand.coef=FALSE){
		if(mode(MAX.p)!="numeric")	stop("MAX.p must be numeric.")
		#if(MAX.p < dfgps_result$p) stop("p must be less than or equal to MAX.p.")
		if(stand.coef!=TRUE && stand.coef!=FALSE) stop('"stand.coef" must be TRUE or FALSE')
		if(intercept!=TRUE && intercept!=FALSE) stop('"intercept" must be TRUE or FALSE')
        #library.dynam("msgps")
                 p = dfgps_result$p
                 p = as.integer(p)
                delta_t=dfgps_result$delta_t
                delta_t=as.double(delta_t)
                betavec = as.integer(dfgps_result$coefficient_index)
                 betaMATRIX=.Call("betaOUT_MATRIX",betavec,p,delta_t)
				 tuning_result <- dfgps_result$tuning
				if(stand.coef==FALSE) {
					a=rep(sqrt(apply(dfgps_result$X,2,var)),ncol(betaMATRIX))
					b=matrix(a,nrow(betaMATRIX),ncol(betaMATRIX))
					betaMATRIX <-  betaMATRIX  / (sqrt(dfgps_result$N-1) * b)
					tuning_result<- dfgps_result$tuning_stand
				}
					if(intercept==TRUE) betaMATRIX <- rbind(rep(mean(dfgps_result$y),ncol(betaMATRIX)),betaMATRIX)
         return(invisible(list(coefficient_mat=betaMATRIX,tuning=tuning_result)))
}

