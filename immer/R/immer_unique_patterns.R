################################################################
# extract unique item response patterns from dataset
immer_unique_patterns <- function( dat , w = rep(1,nrow(dat) ) ){
	I <- ncol(dat)
	N <- nrow(dat)
	v1 <- "P"
	v2 <- NULL
	for (ii in 1:I){
		# ii <- 1
		dig2 <- log( max( dat , na.rm=TRUE ) + 2 , 10 ) >= 1
		na_code <- if ( dig2 ){  99 } else { 9 }
		y <- dat[,ii]
		y[is.na(y)  ] <- na_code
		v1 <- paste0( v1 , y ) 
				}

	f1 <- base::rowsum( matrix(w,ncol=1) , v1 )			
	NP <- nrow(f1)
	y <- matrix( 0 , nrow=NP , ncol=I)
	colnames(y) <- colnames(dat)	
	#*** compute dataframe of indices
	dfr <- data.frame("item" = 1:I )
	h1 <- if( dig2 ){  seq(1,2*I,2) } else { 1:I }
	dfr$start <- h1 + 1
	dfr$end <- dfr$start + dig2 	
	#*** extract data
	for (ii in 1:I){
		y0 <- substring( paste(rownames(f1)) , dfr[ii,"start"] , 
				dfr[ii,"end"] )
		h1 <- as.numeric( y0 )
		h1[ h1 == na_code ] <- NA
		y[,ii ] <- h1
			}
	#*** output			
	res <- list( "y" = y , "w" = f1[,1] , "y_string" = paste( rownames(f1) ) )
	return(res)
		}
################################################################