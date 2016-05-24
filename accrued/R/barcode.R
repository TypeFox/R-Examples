barcode = function(x, ...) { 

	## Throw an error if the argument is not of the correct class.
	if( class(x) != "accrued" )  stop("ERROR: argument is not an object of the 'accrued' class.")

	accrued_data = x		
	data = accrued_data[["data"]]
	STACKED = stackedUploadData( accrued_data )	

	UPLOAD_DATE_vs_NUM_ADDED = as.data.frame( cbind( STACKED[ , "UploadDate"]  , STACKED[ , "NumberAdded"] ) )
	
	
	## X is the vector of upload dates (aka receipt dates).
	X = sort(  unique( STACKED[, "UploadDate"] ) )

	## Y is the vector of upload indicators (2 if upload with positive counts, 1 if "0" counts uploaded, 0 if NA ).
	Y = rep( 0, times = length(X) )

	for( R in 1:(length(X)-1) ) {
		UPLOAD_DATE = X[R]
		temp = UPLOAD_DATE_vs_NUM_ADDED[ UPLOAD_DATE_vs_NUM_ADDED[,1]==UPLOAD_DATE & !is.na(UPLOAD_DATE_vs_NUM_ADDED[,2])  ,  ] 
		if( dim(temp)[1] > 0 ) { 
			temp = temp > 0
			if( sum(temp[,2]) > 0 ) Y[R] = 2 	## At least one count was positive on that upload date.
			else Y[R] = 1	 					## No counts were positive and at least one count zero on that upload date.
		}
	}

	## Set margins.
	BOTTOM_MARGIN = 0.5
	TOP_MARGIN = LEFT_MARGIN = RIGHT_MARGIN = 0.25
 	par( mai=c( BOTTOM_MARGIN, LEFT_MARGIN, TOP_MARGIN, RIGHT_MARGIN ) ) 
	 
	## Set graph limits.
	X_MIN = min(X); X_MAX = max(X); Y_MIN = min(Y);  Y_MAX = max(Y);

	plot( X, Y, 
		  xlim=c(X_MIN,X_MAX), ylim=c(Y_MIN,Y_MAX), 
	      xlab="", ylab="", 
	      main="", type='n',  
	      axes=FALSE,  ...)

	polygon(  c( X_MIN-1, X_MAX+1, X_MAX+1, X_MIN-1 ), 	
			  c( Y_MIN-1, Y_MIN-1, Y_MAX+1, Y_MAX+1 ), 
			  col=rgb(1, .75, .5), border=FALSE ) 
	
	abline( v=(X_MIN-1):(X_MAX+1)+0.5, cex=0.15, col="white" )
	abline( v=X[Y==2], cex=0.50, col="darkblue" )
	abline( v=X[Y==1], cex=0.50, col="gray" )
	

	## x-axis labels.
	JUMP = round( (length(X))/50, 0 ) * 10
	if( JUMP == 0 ) JUMP = 4
	NUMBER_OF_LABELS = round(length(X)/JUMP,0)
	X_TICK_PLACES = 0:NUMBER_OF_LABELS*JUMP
	X_LABELS = X_TICK_PLACES
	LABEL_ORIENTATION = 1

	START_DATE 	= accrued_data[["start"]][[1]]
	if( START_DATE != 1){
		# x-axis labels need to be dates.
		X_LABELS = as.Date(START_DATE + X_TICK_PLACES, origin=as.Date("1970-01-01")) 
		LABEL_ORIENTATION = 2
	}
			
	axis(1, at = X_TICK_PLACES, labels=X_LABELS, cex.axis=0.8, las=LABEL_ORIENTATION, font.axis=1) 
	

} 

