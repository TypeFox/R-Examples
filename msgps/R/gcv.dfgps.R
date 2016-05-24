gcv.dfgps <-
function(dfgps_result,intercept=TRUE, stand.coef=FALSE){
		if(stand.coef!=TRUE && stand.coef!=FALSE) stop('"stand.coef" must be TRUE or FALSE')
		if(intercept!=TRUE && intercept!=FALSE) stop('"intercept" must be TRUE or FALSE')
				  N <- dfgps_result$N
                  gcv_vec <- dfgps_result$RSS / (1-dfgps_result$df / N)^2
                  gcv_vec <- gcv_vec / N
                 gcv_index <- which.min(gcv_vec)
                 gcv_index <- as.integer(gcv_index)
                 p = dfgps_result$p
                 p = as.integer(p)
                 delta_t = dfgps_result$delta_t
                 delta_t=as.double(delta_t)
                betavec = as.integer(dfgps_result$coefficient_index)
				tuning <- dfgps_result$tuning[gcv_index]
  				if(stand.coef==FALSE) tuning<- dfgps_result$tuning_stand[gcv_index]
                 beta_gcv=coef.dfgps(dfgps_result,tuning,intercept,stand.coef)$coefficient
				df <-  dfgps_result$df[gcv_index]
                 return(list(result=gcv_vec,step=gcv_index,tuning=tuning,coef=beta_gcv,df=df))
}

