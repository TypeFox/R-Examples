


#*******************************************************
# summary.btm
summary.btm <- function( object , file=NULL , digits=4,... ){
	
	res <- object

	if ( ! is.null( file ) ){
		sink( paste0( file , "__SUMMARY.Rout") , split=TRUE )
						}	
	
	cat("------------------------------------------------------------\n")
		d1 <- utils::packageDescription("sirt")
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
		cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
		cat("Computation time:" , print(object$s2 - object$s1), "\n\n")

	cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")			
		
    cat("Bradley-Terry Model with Ties and Home Advantage Parameters\n")

	cat("------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
#    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
#    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of individuals = " , object$ic$n , "\n" )    
	cat( "Number of pairwise comparisons = " , object$ic$D , "\n" )    
#    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
	
	cat("------------------------------------------------------------\n")
		cat("Ties and Home advantage parameters\n")
		obji <- res$pars
		V <- ncol(obji)
		for (vv in 3:V ){ 
			obji[,vv] <- round( obji[,vv] , digits ) 
					}
		rownames(obji) <- NULL
		print( obji )   

	cat("------------------------------------------------------------\n")
		cat("Summary of individual effects parameters\n")
		obji <- res$summary.effects
		V <- ncol(obji)
		for (vv in 1:V ){ 
			obji[,vv] <- round( obji[,vv] , digits ) 
					}
		rownames(obji) <- NULL
		print( obji )                

	cat("------------------------------------------------------------\n")
		cat("MLE reliability (separation reliability)\n")
		cat(paste0("MLE Rel = " , round( res$mle.rel , digits )  , "\n") )
		cat(paste0("Separation index = " , round( res$sepG , digits )  , "\n") )

		
	cat("------------------------------------------------------------\n")
		cat("Individual effects parameters\n")
		obji <- res$effects
		V <- ncol(obji)
		for (vv in 3:V ){ 
			obji[,vv] <- round( obji[,vv] , digits ) 
					}
		rownames(obji) <- NULL
		print( obji )                

	# close file
	if ( ! is.null( file ) ){  sink()	}	
		
                }
#*******************************************************



