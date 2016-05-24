
########################################################
# summary tam.modelfit method

summary.tam.modelfit <- function( object , ... ){
  #*****
  cat("Test of Global Model Fit (Maximum Chi Square)\n")
  # cat("p =" , round( object$modelfit.test$p.holm , 5 ) )
  print( round( object$modelfit.test , 5 ) ) 
  
  #******
  cat("\nMADaQ3 Statistic and Test of Global Model Fit (Maximum aQ3)\n")
  obji <- round( object$stat.MADaQ3 , 4)
  print(obji)
  
  #****
  cat("\nSummary of Q3 and adjusted Q3 statistics (based on posterior distribution)\n")
  obji <- object$Q3_summary
  obji[ , -1 ] <- round( obji[,-1] , 4 )
  print(obji)
  
  #*****
  cat("\nFit Statistics\n")
  obji <- object$fitstat
  print( round( obji,3))		
  #*****
  #	cat("\nItem-wise Fit Statistics\n")
  #	obji <- object$chisquare.itemfit
  #	for (vv in 2:(ncol(obji))){ obji[,vv] <- round( obji[,vv] , 3 ) }
  #	print(  obji )	
}
##############################################################

