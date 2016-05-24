uploadPattern =  function(x, horizontal=FALSE, ...) { 
	
	## Throw an error if the argument is not of the correct class.
	if( class(x) != "accrued" )  stop("ERROR: argument is not an object of the 'accrued' class.")

	accrued_data = x		
			
		
	MAX_LAGS = ncol(accrued_data[["data"]]) - 1	
	NCOL = ncol(accrued_data[["data"]])
	NROW = nrow(accrued_data[["data"]])
	STACKED = as.data.frame(stackedUploadData( accrued_data ))

	LAGS = 0:MAX_LAGS

	## X is the vector of upload dates.
	X = STACKED[, "UploadDate"]
	## Y is either vector of encounter dates (horizontal=F) or the vector of lags (horizontal=T).
	Y = STACKED[, "Lag"]

	## This adds horizontal padding on the left-hand side of the plot.
	HORIZONTAL_PADDING=5
	
	## This adds a shift of 1 or zero, depending on the type of plot given.
	LOWER_VERTICAL_SHIFT = 0

	## Specifications if a diagonal plot is desired. 
	## This must come before Y_MIN and Y_MAX are specified since Y is reassigned here.
	ASPECT_RATIO = 1
	if( !horizontal ){
		Y = STACKED[, "EncounterDate"]
		ASPECT_RATIO = 1
		HORIZONTAL_PADDING = 10 	
		LOWER_VERTICAL_SHIFT = 1
	}

	## Plot limits.
	X_MIN = min(X, na.rm=T)
	X_MAX = max(X, na.rm=T)
	Y_MIN = min(Y, na.rm=T) 
	Y_MAX = max(Y, na.rm=T)

	if( X_MIN == X_MAX ) X_MAX = X_MIN + 1
	if( Y_MIN == Y_MAX ) Y_MAX = Y_MIN + 1

	## Aspect Ratio is optimized for around 100-200 upload dates and 20 lags.
	if( horizontal ) { ASPECT_RATIO = (Y_MAX-Y_MIN)/(X_MAX-X_MIN) + 0.025 }
	X_LABEL = Y_LABEL = ""
	BOTTOM_MARGIN = LEFT_MARGIN = 0.5
	TOP_MARGIN = RIGHT_MARGIN = 0.25

	## If a start date is provided, make margins large enough to accommodate the date labels.
	START_DATE = accrued_data[["start"]][[1]]
	if( START_DATE != 1  )  BOTTOM_MARGIN = LEFT_MARGIN = 1
	
 	par(  mai=c( BOTTOM_MARGIN, LEFT_MARGIN, TOP_MARGIN, RIGHT_MARGIN )  )
 	plot( X, Y, 
 	      xlim=c(X_MIN, X_MAX), 
 	      ylim=c(Y_MIN, Y_MAX+1), 
 	      xlab=X_LABEL, 
 	      ylab=Y_LABEL,
	      main="", type='n', 
	      axes=FALSE, 
	      xaxs="i", yaxs="i" , ...)

	## Draw background orange rectangle.
	if( horizontal )	{
		polygon(  c( X_MIN-HORIZONTAL_PADDING  , X_MAX+MAX_LAGS			   , X_MAX+MAX_LAGS, X_MIN-HORIZONTAL_PADDING ), 	
			      c( Y_MIN-LOWER_VERTICAL_SHIFT, Y_MIN-LOWER_VERTICAL_SHIFT, Y_MAX+1       , Y_MAX+1 ), 
			      col=rgb(1, 0.75, 0.5), border=FALSE ) 
	}
	else {
		polygon(  c( X_MIN-HORIZONTAL_PADDING  , X_MAX+HORIZONTAL_PADDING  , X_MAX+HORIZONTAL_PADDING, X_MIN-HORIZONTAL_PADDING ), 	
			      c( Y_MIN-LOWER_VERTICAL_SHIFT, Y_MIN-LOWER_VERTICAL_SHIFT, Y_MAX+1                 , Y_MAX+1 ), 
			      col=rgb(1, 0.75, 0.5), border=FALSE ) 
	}


	####################################################################
	## Draw the dark blue box at a (i,j) if an upload on day i 
	## containing data from day j was uploaded, with positive counts.
	## Draw the gray box at a (i,j) if an upload on day i containing 
	## data from day j was uploaded, with zero counts.
	## Otherwise no box is drawn.
	####################################################################

	## This is the vector of numbers added.
	NUM_ADDED = STACKED[ , "NumberAdded"]
	NOT_NA_AND_POS  = (NUM_ADDED > 0)  & !is.na(NUM_ADDED)
	NOT_NA_AND_ZERO = (NUM_ADDED == 0) & !is.na(NUM_ADDED)

	NO_PRIOR_DATA = rep(T, times=length(NUM_ADDED)) 
	for( R in 1:length(NO_PRIOR_DATA) ) {
		lag = STACKED[R, "Lag"]
		encounterDate = STACKED[R, "EncounterDate"]
		if( lag > 0 ){ 
			temp = STACKED[ STACKED$EncounterDate==encounterDate & STACKED$Lag<lag & !is.na(STACKED$Counts),"Counts"] 
			if( length(temp) > 0 )  
				if( sum(temp > 0) > 0) NO_PRIOR_DATA[R] = F
		}
		
	}
		
	## Positive counts are set to 2, "0" counts, or <0 counts remain at 1.
	## But, if counts have already been received for that encounter date and lag combination,
	SELECTION_VECTOR = rep(0, times=length(NUM_ADDED)) 
	for( R in 1:length(SELECTION_VECTOR) ) {
		if( NOT_NA_AND_POS[R] ) SELECTION_VECTOR[R] = 2
		if( NOT_NA_AND_ZERO[R]& NO_PRIOR_DATA[R] ) SELECTION_VECTOR[R] = 1
	}

	for( R in 1:length(SELECTION_VECTOR)) { 
		if( SELECTION_VECTOR[R] == 2 ) polygon( c( X[R], X[R]+1, X[R]+1, X[R]   ), 
                                            c( Y[R], Y[R]  , Y[R]+1, Y[R]+1 ), 
                                            col="darkblue", border=FALSE ) 
		if( SELECTION_VECTOR[R] == 1 ) polygon( c( X[R], X[R]+1, X[R]+1, X[R]   ), 
                                            c( Y[R], Y[R]  , Y[R]+1, Y[R]+1 ), 
                                            col="gray", border=FALSE ) 
 	}
			
	## Add vertical white lines to make the boxes more distinct.
	if( horizontal ) abline( v=(X_MIN-HORIZONTAL_PADDING):(X_MAX+MAX_LAGS), cex=0.15, col="white" ) 	
	else abline( v=(X_MIN-HORIZONTAL_PADDING):(X_MAX+HORIZONTAL_PADDING), cex=0.15, col="white" ) 	
	
	## Add horizontal white lines to make the boxes more distinct.
	abline( h=(Y_MIN-1):(Y_MAX+1), cex=0.15, col="white" )		

	## Add a black horizontal line of height 0.
	abline( h=0, cex=0.30, col="black" )		


	## Add labels.	
	X_NUMBER_OF_LABELS = 10
	X_JUMP = round( (NROW/X_NUMBER_OF_LABELS)/10, 0 )*10
	if( X_JUMP == 0 ) X_JUMP = 4
	X_TICK_PLACES	 = (0:(X_NUMBER_OF_LABELS)) * X_JUMP
	X_LABELS	 =  X_TICK_PLACES
	X_LABEL_ORIENTATION = 1

	Y_JUMP = round(MAX_LAGS/4, 0)
	if( Y_JUMP == 0 ) Y_JUMP = 4
	## These are the plot settings if "horizontal==T" and there are more than 10 lags.
	Y_TICK_PLACES = seq(0,MAX_LAGS, by=Y_JUMP)
	Y_LABELS = Y_TICK_PLACES
	Y_CEX_AXIS = 0.8
	Y_FONT_AXIS = 1


	## If a start date is specified, set the x-labels which are the same regardless of the value of "horizontal".
	if( START_DATE != 1 ) {
		# x-axis labels need to be dates.
		# x-axis labels are the same, regardless of whether "horizontal" is TRUE or FALSE.
		X_LABELS = as.Date(START_DATE + X_TICK_PLACES, origin=as.Date("1970-01-01")) 
		X_LABEL_ORIENTATION = 2
	} 

	## The y-labels differ depending on whether "horiztonal" is T or F. If it is T, then either
	## the default settings from above or used, or they're tweaked in the event of a low value of MAX_LAGS.
	## Otherwise (if "horizontal==F"), the y-labels are the same as the x-labels.
	if( horizontal ) { 
		if( Y_MAX-Y_MIN <= 10 ) { 
			Y_TICK_PLACES = Y_LABELS = c(0,5,10)
			Y_CEX_AXIS = 0.4
			axis(2, at=Y_TICK_PLACES, labels=Y_LABELS, las=2, cex.axis= Y_CEX_AXIS, font.axis=Y_FONT_AXIS)
		}
	} else {
		Y_TICK_PLACES = X_TICK_PLACES
		Y_LABELS = X_LABELS
		abline(h=0, col="darkgrey", cex=1)
	}

	axis(1, at=X_TICK_PLACES, labels=X_LABELS, las=X_LABEL_ORIENTATION, cex.axis=0.8,        font.axis=1) 
	axis(2, at=Y_TICK_PLACES, labels=Y_LABELS, las=2,                 , cex.axis=Y_CEX_AXIS, font.axis=Y_FONT_AXIS)


} 
	
