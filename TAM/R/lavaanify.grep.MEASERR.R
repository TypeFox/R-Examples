
#####################################################
# grep term "MEASERR" for measurement errors
lavaanify.grep.MEASERR <- function( lavmodel ){

	lavmodel0 <- lavmodel
	lavmodel <- gsub( ";" , "\n" , lavmodel , fixed=TRUE)
	lavmodel <- gsub( " " , "" , lavmodel )
	syn <- strsplit( lavmodel , split="\n" , fixed=TRUE )[[1]]
	syn <- syn[ syn != "" ]
	SS <- length(syn)
	dfr <- data.frame( "index" = 1:SS , "syn" = syn )
	dfr$MEASERR <- 1 * ( substring( dfr$syn , 1 , 7 ) == "MEASERR" )
	dfr$MEASERR1 <- 1 * ( substring( dfr$syn , 1 , 8 ) == "MEASERR1" )	
	dfr$MEASERR0 <- 1 * ( substring( dfr$syn , 1 , 8 ) == "MEASERR0" )		
	dfr$true <- NA
	dfr$obs <- NA
	dfr$errvar <- NA

	ind <- which( dfr$MEASERR == 1 )
	if ( length(ind) > 0 ){ 
	# decompose term MEASERR
	  for (ii in ind ){
		#ii <- ind[1]
		l1 <- paste(dfr$syn[ii] )
		l1 <- gsub( "MEASERR1(" , "" , l1 , fixed=TRUE )
		l1 <- gsub( "MEASERR0(" , "" , l1 , fixed=TRUE )
		l1 <- gsub( ")" , "" , l1 , fixed=TRUE )
		l2 <- strsplit( l1 , "," )[[1]]
		dfr[ ii , "true" ] <- l2[1]
		dfr[ ii , "obs" ] <- l2[2]
		dfr[ ii , "errvar" ] <- as.numeric( l2[3] )
						}
	 # create new string of lavaan syntax
	  N1 <- nrow(dfr)
	  lavmodel0 <- "" 
	  for (nn in 1:N1){
		if ( dfr$MEASERR[nn] == 0 ){
            lavmodel0 <- paste0( lavmodel0 , "\n" , dfr[nn,"syn"] )		
								}
		if ( dfr$MEASERR[nn] == 1 ){
            lavmodel0 <- paste0( lavmodel0 , "\n" , 
							dfr[nn,"true"] , "=~1*" , dfr[nn,"obs"] )		
            lavmodel0 <- paste0( lavmodel0 , "\n" , 
							dfr[nn,"obs"] , "~~" , dfr[nn,"errvar"] , "*" , dfr[nn,"obs"] )									
								}												
		if ( dfr$MEASERR0[nn] == 1 ){
            lavmodel0 <- paste0( lavmodel0 , "\n" , 
							dfr[nn,"true"] , "~~" , dfr[nn,"true"]  )									
								}												
								
					}	
					}
	return(lavmodel0)
		}
########################################################################		
