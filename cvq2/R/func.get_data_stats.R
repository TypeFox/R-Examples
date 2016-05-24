func.get_data_stats <-
function( data, tmp, stat_type, nu ){
  dt_stat <- NULL
  
  # mean values
  dt_stat$observed_mean <- mean( data[,tmp$obs_col] )
  dt_stat$predicted_mean <- mean( data[,tmp$pred_col] )
  
  observed_mean = dt_stat$observed_mean 
  
  if( stat_type == "prediction" )
    observed_mean = data[,tmp$obs_mean_col]
    
  dt_stat$x2 <- func.calcXSquare( data[,tmp$pred_col], data[,tmp$obs_col], observed_mean )

  # rmse
  #   nu: degrees of freedom (in general)
  #   linear regression, prediction with external data set == 0
  #   cross validation == 1 (control/random sample - Stichprobe)
  dt_stat$rmse <- func.calcRMSE( data[,tmp$pred_col], data[,tmp$obs_col], nu )
  
  return( dt_stat )
}

