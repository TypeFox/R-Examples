##########################################################################################

plot.accrued =  function(x, ...) { 

	## Throw an error if the argument is not of the correct class.
	if( class(x) != "accrued" )  stop("ERROR: argument is not an object of the 'accrued' class.")

	accrued_data = x		

	## Get largest lag from data (lag is counted from zero, so it will be maximum column number minus 1).
	MAX_LAGS = ncol( accrued_data[["data"]] ) - 1

	## Throw an error if there are not enough columns of data.
	if( MAX_LAGS == 0 ) stop("The function 'plot.accrued' is useful only if there are at least two lag columns.")

	NROW = nrow( accrued_data[["data"]] )

	## The "stackedUploadData" takes the "data" component of a data.accrued object and
	## stacks it, adding "lag" as a column variable. It also calculates the number of counts added
	## on a particular encounter date from one upload date to the next.
	STACKED = stackedUploadData( accrued_data )
 	
 
 	###############################################################################################
 	###############################################################################################
 	###############################################################################################
 	
	LAYER_LABELS = c("0")
	for( lag in 1:MAX_LAGS ) LAYER_LABELS = rbind( LAYER_LABELS ,paste( lag-1, "-", lag, sep="" ) )
	LAYER_LABELS = as.vector(LAYER_LABELS)
	COL_LABELS = c( "minVal", 
					"maxVal", 
					"diff",
					"grayPadding", 
					"whitePadding", 
					"relLow", 
					"relZero", 
					"relHigh", 
					"absGrayBoxLow", 
					"absGrayBoxHigh",
					"absGraphLow",
					"absGraphZero",
					"absGraphHigh" ) 
 	LAYER_INFO = matrix(  0, nrow=length(LAYER_LABELS), ncol=length(COL_LABELS), 
					   dimnames=list( LAYER_LABELS, COL_LABELS) )	
	LAYER_INFO = as.data.frame(LAYER_INFO)
 	for( lag in 0:MAX_LAGS ) {
 		LAYER_INFO[lag+1,"minVal"] = min(STACKED$NumberAdded[STACKED$Lag == lag], na.rm=T)
 		LAYER_INFO[lag+1,"maxVal"] = max(STACKED$NumberAdded[STACKED$Lag == lag], na.rm=T) 		
 		LAYER_INFO[lag+1,"diff"] = LAYER_INFO[lag+1,"maxVal"]  - LAYER_INFO[lag+1,"minVal"]
 	}
	temp = LAYER_INFO[ !is.na(LAYER_INFO$diff), "diff"]
	temp = temp[temp>0]
	padding = 5
	if(length(temp)>0) padding = max(temp)/4  

### STOPPED HERE	
### STOPPED HERE	
### STOPPED HERE	

	LAYER_INFO[ ,"grayPadding"] = rep(padding/2, times=length(LAYER_LABELS))
	LAYER_INFO[ ,"whitePadding"] = rep(padding, times=length(LAYER_LABELS))

	###################################################
	## this is a temporary proof-of-concept for now. ##
	###################################################	
	for( lag in 0:MAX_LAGS ) {
		if(LAYER_INFO[lag+1,"minVal"] == LAYER_INFO[lag+1,"maxVal"])	 {
			LAYER_INFO[lag+1,"relHigh"] = LAYER_INFO[lag+1,"grayPadding"] 
			LAYER_INFO[lag+1,"relLow"]  = -LAYER_INFO[lag+1,"relHigh"] 
		} else if(LAYER_INFO[lag+1,"minVal"] >= 0 & LAYER_INFO[lag+1,"maxVal"] >= 0 )	 {
			LAYER_INFO[lag+1,"relHigh"] = LAYER_INFO[lag+1,"grayPadding"] 
			LAYER_INFO[lag+1,"relLow"]  = -LAYER_INFO[lag+1,"relHigh"] 
		} else if(LAYER_INFO[lag+1,"minVal"] <= 0 & LAYER_INFO[lag+1,"maxVal"] <= 0 )	 {
			LAYER_INFO[lag+1,"relHigh"] = LAYER_INFO[lag+1,"grayPadding"] 
			LAYER_INFO[lag+1,"relLow"]  = -LAYER_INFO[lag+1,"relHigh"] 
		} else {
			LAYER_INFO[lag+1,"relHigh"] = LAYER_INFO[lag+1,"maxVal"] 
			LAYER_INFO[lag+1,"relLow"]  = LAYER_INFO[lag+1,"minVal"] 
		}
	}
	##########
	lag = 0
	LAYER_INFO[lag+1,"absGrayBoxLow"] =  0 + LAYER_INFO[lag+1,"whitePadding"]
	LAYER_INFO[lag+1,"absGraphZero"] = ( LAYER_INFO[lag+1,"absGrayBoxLow"] +	LAYER_INFO[lag+1,"grayPadding"]
										+	abs(LAYER_INFO[lag+1,"minVal"])   )
	LAYER_INFO[lag+1,"absGrayBoxHigh"] = ( LAYER_INFO[lag+1,"absGraphZero"] + LAYER_INFO[lag+1,"maxVal"]
										+	abs(LAYER_INFO[lag+1,"grayPadding"])	)
	LAYER_INFO[lag+1,"absGraphLow"] = ( LAYER_INFO[lag+1,"absGraphZero"] + LAYER_INFO[lag+1,"relLow"] )
	LAYER_INFO[lag+1,"absGraphHigh"] = ( LAYER_INFO[lag+1,"absGraphZero"] + LAYER_INFO[lag+1,"relHigh"] )
	for( lag in 1:MAX_LAGS ) {
		LAYER_INFO[lag+1,"absGrayBoxLow"] =  LAYER_INFO[lag,"absGrayBoxHigh"] + LAYER_INFO[lag+1,"whitePadding"]
		LAYER_INFO[lag+1,"absGraphZero"] = ( 	
											LAYER_INFO[lag+1,"absGrayBoxLow"]
										+	LAYER_INFO[lag+1,"grayPadding"]
										+	abs(LAYER_INFO[lag+1,"minVal"])
										)
		LAYER_INFO[lag+1,"absGrayBoxHigh"] = ( 	
											LAYER_INFO[lag+1,"absGraphZero"]
										+	LAYER_INFO[lag+1,"maxVal"]
										+	abs(LAYER_INFO[lag+1,"grayPadding"])
										)
		LAYER_INFO[lag+1,"absGraphLow"] = ( LAYER_INFO[lag+1,"absGraphZero"] + LAYER_INFO[lag+1,"relLow"] )
		LAYER_INFO[lag+1,"absGraphHigh"] = ( LAYER_INFO[lag+1,"absGraphZero"] + LAYER_INFO[lag+1,"relHigh"] )
	}
	##############################################################################################################

 	FIRST_NON_NA_ENCOUNTER_DATE = min( STACKED[ !is.na(STACKED[,"NumberAdded"]) , "EncounterDate" ] )
 	FINAL_NON_NA_ENCOUNTER_DATE = max( STACKED[ !is.na(STACKED[,"NumberAdded"]) , "EncounterDate" ] )
	X = FIRST_NON_NA_ENCOUNTER_DATE: FINAL_NON_NA_ENCOUNTER_DATE
	if( length(X) <= 1 ) stop("Not encounter dates with non-missing data.")
	Y = rep(0, times=length(X))
	Y[length(Y)] = LAYER_INFO[MAX_LAGS+1,"absGrayBoxHigh"] 
	X_MIN = min(X)
	X_MAX = max(X)
	Y_MIN = min(Y)
	Y_MAX = max(Y)
	plot( X, Y, 
		  xlim=c(X_MIN,X_MAX), ylim=c(Y_MIN,Y_MAX), 
	      xlab="", ylab="", 
	      main="", type='n',  
	      axes=FALSE, xaxs="i", yaxs="i" )#,  ...)
	for( lag in 0:MAX_LAGS ) { 
		polygon(  c( X_MIN-1, X_MAX+1, X_MAX+1, X_MIN-1 ), 	
				  c( LAYER_INFO[lag+1,"absGrayBoxLow"], LAYER_INFO[lag+1,"absGrayBoxLow"], 
				  		LAYER_INFO[lag+1,"absGrayBoxHigh"], LAYER_INFO[lag+1,"absGrayBoxHigh"] ), 
				  col=rgb(0.9, 0.9, 0.9), border=FALSE ) 
		# rgb(1, .95, .85)
		abline( h=LAYER_INFO[lag+1,"absGraphZero"] ,col=rgb(0.3, 0.3, 0.3) ) # cex=0.15, 
	}
	#########################################################################################################################
 	### y-axis labels.
 	LAYER_LABELS
 	Y_AXIS_INDEX 		= (0: MAX_LAGS)+1
	## Y_AXIS_INDEX  		= seq( 0, MAX_LAGS, by=min(4, floor(MAX_LAGS/2)) ) + 1
	Y_LABEL_TEXT   		= LAYER_LABELS[Y_AXIS_INDEX] 
	Y_LABEL_HEIGHTS  	= LAYER_INFO[Y_AXIS_INDEX,"absGraphZero"] 
	axis( 2, at=Y_LABEL_HEIGHTS, labels=Y_LABEL_TEXT, cex.axis=0.8, las=2, font.axis=1 )

 	### x-axis labels.
 	NUMBER_OF_LABELS = 10
	JUMP = round( (NROW/NUMBER_OF_LABELS)/5, 0 ) * 10
	if( JUMP == 0 ) JUMP = 2
	X_TICK_PLACES = ( 0:NUMBER_OF_LABELS ) * JUMP
	X_LABELS = X_TICK_PLACES
	LABEL_ORIENTATION = 1
	
	START_DATE 	= accrued_data[["start"]][[1]]
	if( START_DATE != 1){
		# x-axis labels need to be dates.
		X_LABELS = as.Date(START_DATE + X_TICK_PLACES, origin=as.Date("1970-01-01")) 
		LABEL_ORIENTATION = 2
	}

	abline(v=1,col="black", lwd=2)			
	axis(1, at = X_TICK_PLACES, labels=X_LABELS, cex.axis=0.8, las=LABEL_ORIENTATION, font.axis=1) 

	#########################################################################################################################

	X_VALUES = X
	lag = 0 
	COLOR="blue"
	rel_Y_VALUES = STACKED[ (STACKED[,"Lag"] == lag) , "NumberAdded" ]	 
	Y_VALUES = rel_Y_VALUES + LAYER_INFO[lag+1,"absGraphZero"]
	points( X_VALUES, Y_VALUES, pch=".", cex=0.5, col=COLOR )
	lines(  X_VALUES, Y_VALUES, col=COLOR)

	for( lag in 1:(MAX_LAGS) ) {
		COLOR="blue"
		if( lag %% 2 == 0 ) COLOR = "forestgreen"
		rel_Y_VALUES = STACKED[ (STACKED[,"Lag"] == lag) , "NumberAdded" ]	 
		Y_VALUES = rel_Y_VALUES + LAYER_INFO[lag+1,"absGraphZero"]
		points( X_VALUES, Y_VALUES, pch=".", cex=0.5, col=COLOR )
		lines(  X_VALUES, Y_VALUES, col=COLOR, lwd=1)
	}

} 



