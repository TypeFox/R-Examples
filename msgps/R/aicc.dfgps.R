aicc.dfgps <-
function(dfgps_result,intercept=TRUE, stand.coef=FALSE){
		if(stand.coef!=TRUE && stand.coef!=FALSE) stop('"stand.coef" must be TRUE or FALSE')
		if(intercept!=TRUE && intercept!=FALSE) stop('"intercept" must be TRUE or FALSE')
                  N <- dfgps_result$N
				tau2 <-  dfgps_result$RSS / N
                  aicc_vec <- N*log(2*pi*tau2) + N + 2*N*dfgps_result$df / (N-dfgps_result$df-1)
                 aicc_index <- which.min(aicc_vec)
                 aicc_index <- as.integer(aicc_index)
                 p = dfgps_result$p
                 p = as.integer(p)
                 delta_t = dfgps_result$delta_t
                 delta_t=as.double(delta_t)
                betavec = as.integer(dfgps_result$coefficient_index)
				tuning <- dfgps_result$tuning[aicc_index]
  				if(stand.coef==FALSE) tuning<- dfgps_result$tuning_stand[aicc_index]
               beta_aicc=coef.dfgps(dfgps_result,tuning,intercept,stand.coef)$coefficient
				df <-  dfgps_result$df[aicc_index]
                 return(list(result=aicc_vec,step=aicc_index,tuning=tuning,coef=beta_aicc,df=df))
}

