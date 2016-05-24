cp.dfgps <-
function(dfgps_result,tau2,intercept=TRUE, stand.coef=FALSE){
		if(stand.coef!=TRUE && stand.coef!=FALSE) stop('"stand.coef" must be TRUE or FALSE')
		if(intercept!=TRUE && intercept!=FALSE) stop('"intercept" must be TRUE or FALSE')
				X <- dfgps_result$X
				X <- scale(X) / sqrt(nrow(X)-1)	
				X[is.nan(X)] <- 0
				y <- dfgps_result$y
				y <- y-mean(y)
				if(missing(tau2)){
				if( nrow(X) <= ncol(X)){
						tau2 <- var(y)
						}else if(det(t(X)%*%X) < 1e-3){
						tau2 <- sum((y-X%*%solve(t(X)%*%X + 0.001*diag(ncol(X))) %*%t(X)%*%y)^2) / (length(y) - ncol(X))
					}else{
						tau2 <- sum((y-X%*%solve(t(X)%*%X) %*%t(X)%*%y)^2) / (length(y) - ncol(X))
					}
				}
				if(mode(tau2)!="numeric")	stop("tau2 must be numeric.")
				if(tau2 < 0)	stop("tau2 must be positive.")
                  Mallows_vec <- dfgps_result$RSS / tau2 + 2*dfgps_result$df
                 Mallow_index <- which.min(Mallows_vec)
                 Mallow_index <- as.integer(Mallow_index)
                 p = dfgps_result$p
                 p = as.integer(p)
                 delta_t = dfgps_result$delta_t
                 delta_t=as.double(delta_t)
                betavec = as.integer(dfgps_result$coefficient_index)
				tuning <- dfgps_result$tuning[Mallow_index]
				if(stand.coef==FALSE) tuning<- dfgps_result$tuning_stand[Mallow_index]
				 beta_Mallows=coef.dfgps(dfgps_result,tuning,intercept,stand.coef)$coefficient
				df <-  dfgps_result$df[Mallow_index]
                 return(list(result=Mallows_vec,step=Mallow_index,tuning=tuning,coef=beta_Mallows,df=df))
}

