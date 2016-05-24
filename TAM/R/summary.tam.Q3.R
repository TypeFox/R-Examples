
########################################################
# summary tam.modelfit method

summary.tam.Q3 <- function( object , ... ){
  #****
  cat("Summary of Q3 and adjusted Q3 statistics (based on WLEs)\n")
  obji <- object$Q3_summary
  obji[ , -1 ] <- round( obji[,-1] , 4 )
  print(obji)

}
##############################################################

