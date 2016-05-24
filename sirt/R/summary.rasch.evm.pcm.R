#*******************************************************
# Summary for rasch.mirtlc object                         *
summary.rasch.evm.pcm <- function( object,... ){
    # object      ... object from rasch.mml                #
    cat("------------------------------------------------------------------------------------\n")
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
#	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
#	cat("Computation time:" , print(object$s2 - object$s1), "\n\n")
    cat("Partial Credit Model (estimated with eigenvector method) \n\n")
    options(scipen=999)
	cat( "Number of groups:" , object$desc$G[1] , "\n" )	
	cat( "Number of persons:" , object$desc$Nobs , "\n" )	
	cat( "Weighted number of persons:" , round( (object$desc$sum.weights) , 2 ) , "\n" )				
	cat( "Number of items:" , object$desc$N.items[1] , "\n" )	
	cat( paste0("Number of items per parameter: " , 
			"M = " , round( object$desc$M.Nitems[1],1) , 
			" | SD = " , round( object$desc$SD.Nitems[1],1) , 			
			" | Min = " , object$desc$min.Nitems[1] , 
			" | Max = " , object$desc$max.Nitems[1] , 			
			"\n" ) 
				)		
	cat( "Number of parameters:" , sum(object$desc$N.pars) , "\n\n" )		

	
    cat("------------------------------------------------------------------------------------\n")
	cat("Item Parameters \n\n")
	obji <- object$item
	for (vv in seq(4,ncol(obji)) ){
		obji[,vv] <- round( obji[,vv] , 3 )
					}
	print( obji )      
    if (! is.null( object$difstats) ){	
		cat("\n------------------------------------------------------------------------------------\n")
		cat("DIF Tests \n\n")
		obji <- object$difstats
		for (vv in seq(2,ncol(obji)) ){
			obji[,vv] <- round( obji[,vv] , 3 )
						}
		print( obji )      
			}	
    options(scipen=0)	
			}
#*******************************************************
