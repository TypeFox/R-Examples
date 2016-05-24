

immer_reshape_wideformat <- function( y , pid , rater , Nmin_ratings = 1 ){	
	y_dfr <- FALSE
	if ( is.data.frame(y) ){
		y_dfr <- TRUE			
						}
	if ( ! y_dfr ){
		dfr1 <- immer_reshape_wideformat_vector(y, pid,rater, Nmin_ratings )
				}
	if ( y_dfr ){
		NV <- ncol(y)		
		for (vv in 1:NV){		
			y1 <- as.vector(y[,vv])
			dfr2 <- immer_reshape_wideformat_vector(y=y1,  pid=pid,
					rater=rater, Nmin_ratings=Nmin_ratings )		
			colnames(dfr2)[-1] <- paste0( colnames(y)[vv] , "_" , colnames(dfr2)[-1] ) 
			if ( vv == 1 ){ dfr1 <- dfr2 }
			if ( vv > 1 ){
				dfr1 <- base::merge( x = dfr1 , y = dfr2 , by = "pid" , all=TRUE )
							}
						}
			# dfr1 <- dfr1[ order( dfr1$pid ) , ]			
				}	
	return(dfr1)
		}
		

immer_reshape_wideformat_vector <- function(y, pid,rater, Nmin_ratings ){		
	rater <- paste(rater)
	Nobs <- rowsum( 1 - is.na(y) ,  pid  )
	Nobs <- Nobs[ Nobs >= Nmin_ratings , ]
	persons <- names(Nobs)
	NP <- length(persons)
	data <- data.frame( "pid" = pid , "rater" = rater , "y" = y )
	data <- data[ data$pid %in% persons , ]
	raters <- sort( unique( paste(data$rater )))
	RR <- length(raters)
	y <- matrix( NA , nrow=NP , ncol=RR+1 )
	y <- as.data.frame(y)
	colnames(y) <- c("pid" , raters )
	y$pid <- persons
	indM <- cbind( match( paste(data$pid) , persons) , 
				match( paste(data$rater) , raters)+1 ) 
	y[ indM ] <- data$y
	return(y)		
			}