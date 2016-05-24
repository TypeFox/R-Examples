
#############################################
# summary lsem
summary.lsem <- function( object , file=NULL , digits=3 , ... ){

	# open sink for a file
	osink( file=file , prefix="__SUMMARY.Rout" )

	cat("-----------------------------------------------------------------\n")
	cat("Local Structural Equation Model \n\n")
		
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , " | " )	
    d1 <- utils::packageDescription("lavaan")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )	
	
	
	cat(paste0("Function 'lsem.estimate', type='" , object$type,"'") , "\n\n")
	
	
	cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")		
	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
    
	
	cat( paste0( "Number of observations = " , round(object$N,digits) ) , "\n")
	if ( object$type == "LSEM"){
		cat( paste0( "Bandwidth factor = " , round(object$h,digits) ) , "\n")
		cat( paste0( "Bandwidth = " , round(object$bw,digits) ) , "\n")
		cat( paste0( "Number of focal points for moderator = " , 
							length(object$moderator.grid ) ) , "\n")
								}

	if ( object$type == "MGM"){
#		cat( paste0( "Bandwidth factor = " , round(object$h,digits) ) , "\n")
#		cat( paste0( "Bandwidth = " , round(object$bw,digits) ) , "\n")
		cat( paste0( "Number of groups for moderator = " , 
							length(object$moderator.grid ) ) , "\n")
								}								
								
	cat("\nlavaan Model\n")
	cat(object$lavmodel)						
						
    cat("\n\n")
	cat("Parameter Estimate Summary\n\n")
	obji <- object$parameters_summary
	VV <- ncol(obji)
	for (vv in 2:VV){
		obji[,vv] <- round( obji[,vv] , digits )
					}
	print(obji)
	
    cat("\n")
	cat("Distribution of Moderator: Density and Effective Sample Size\n\n")
	
	obji <- object$moderator.density
	VV <- ncol(obji)
	for (vv in 1:VV){
		obji[,vv] <- round( obji[,vv] , digits )
					}
	print(obji)		
	cat("\n")
	obji <- object$moderator.stat	
	VV <- ncol(obji)
	for (vv in 2:VV){
		obji[,vv] <- round( obji[,vv] , digits )
					}
	print(obji)	
	
	# close file
	csink(file)
	
	}
#############################################