
##############################################################
summary.BIFIEdata <- function(object , ... ){
	cat("------------------------------------------------------------\n")
	d1 <- utils::packageDescription("BIFIEsurvey")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )	
	cat( paste0("Function '" , class(object) ) )
	cat("'" )		
	if( object$cdata ){ cat("\nCompact BIFIEdata") }
	cat("\n\n" )	
	cat( "Date of Analysis:" , paste( object$time[1] ) , "\n\n" )
	
	if ( ! object$NMI ){
		cat("Multiply imputed dataset\n\n")
					}
	if ( object$NMI ){
		cat("Nested multiply imputed dataset\n\n")
		# object$Nimp <- object$Nimp_NMI
					}
					
    cat( "Number of persons =" , object$N , "\n" )    
    cat( "Number of variables =" , object$Nvars , "\n" )
	if ( ! object$NMI){ 	
		cat( "Number of imputed datasets =" , object$Nimp , "\n" ) 
					}
	if ( object$NMI){ 	
		cat( "Number of imputed between-nest datasets =" , object$Nimp_NMI[1] , "\n" ) 
		cat( "Number of imputed within-nest datasets =" , object$Nimp_NMI[2] , "\n" ) 
					}
					
    cat( "Number of Jackknife zones per dataset =" , object$RR , "\n" ) 
	
	fayfac <- object$fayfac
	if ( length(fayfac) == 1){	
		cat( "Fay factor =" , round( object$fayfac , 5 ) , "\n\n" ) 	
							} else {
		mf <- mean(fayfac)
		sdf <- stats::sd(fayfac)		
		cat( "Fay factor: M =" , round( mf , 5 ) , "| SD =" , 
			round( sdf , 5 ) , "\n\n" ) 													
							}
	
	x2 <- BIFIE_object_size(object)
	cat( "Object size: \n  " ) 			
	cat( "Total: " , x2$value_string , "\n")
	
	obj1 <- "datalistM"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")
    
	obj1 <- "datalistM_ind"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")

	obj1 <- "datalistM_imputed"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")
	
	obj1 <- "datalistM_impindex"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")
	
	obj1 <- "dat1"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")	

	obj1 <- "wgtrep"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")	

		
						}
######################################################################	

##################################################
# output object sizes in MB and GB
BIFIE_object_size <- function( x1 ){
	x1 <- utils::object.size(x1)
	vals <- c( x1 , as.numeric(x1 / 1024^2) ,
					as.numeric(x1 / 1024^3)  )
	names(vals) <- c("B" , "MB" , "GB" )		
    res <- list( "value" = vals )
	res$value_string <- paste0( round( vals[2] , 3 ) , 
			" MB | " , round( vals[3] , 5 ) ,	" GB" ) 		
    return(res)
			}