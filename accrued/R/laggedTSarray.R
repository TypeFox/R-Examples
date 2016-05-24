############################################################################## 
## LAGGED TIME SERIES ARRAY.
## This script is for plotting running averages and errors.  
##############################################################################
 	
laggedTSarray =  function(x, daysOfHistory = NULL, lags = NULL, ...) { 
	
	## Throw an error if the argument is not of the correct class.
	if( class(x) != "accrued" )  stop("ERROR: argument is not an object of the 'accrued' class.")

	accrued_data = x
	data = accrued_data[["data"]]
	MAX_LAGS = ncol(data) - 1
	NROW = nrow(data)
	NCOL = ncol(data)

	## HISTORY_DAYS is the number of prior days used in the running median calculations.
	HISTORY_DAYS = NULL
	if( is.null(daysOfHistory) ) HISTORY_DAYS = 30
	else HISTORY_DAYS =  daysOfHistory
	
	if( HISTORY_DAYS < 1 )  stop("ERROR: daysOfHistory must be an integer at least equal to 1.")

	LAGS = lags
	if( is.null(lags) ) LAGS = 0:MAX_LAGS

	# There are three data structures that have to be created:
	# 1. RUNNING_MEDIANS
	# 2. ABSOLUTE_DEVIATIONS
	# 3. MEDIAN_ABSOULTE_DEVIATIONS

 	#####################
	## RUNNING MEDIANS ##
	#####################
	RUNNING_MEDIANS = matrix(NA, ncol=NCOL, nrow=NROW, dimnames=dimnames(data))
 	for( L in LAGS ) 
		for( R in (HISTORY_DAYS+1):NROW ) {
			COL_L = paste(L, sep="")
			TEMP = data[(R-HISTORY_DAYS):(R-1), COL_L]
			if(  sum( !is.na( TEMP ) ) >= HISTORY_DAYS/2  ) 
				RUNNING_MEDIANS[R, COL_L] = median( TEMP, na.rm=T )
		}

	#########################
	## ABSOLUTE DEVIATIONS ##
	#########################
	ABSOLUTE_DEVIATIONS = matrix( NA, ncol=NCOL, nrow=NROW, dimnames=dimnames(data))
 	for( L in LAGS ) 
		for( R in (HISTORY_DAYS+1):NROW ) {
			COL_L = paste(L, sep="")
			if( !is.na( RUNNING_MEDIANS[R, COL_L] ) ) 
				ABSOLUTE_DEVIATIONS[R, COL_L] = abs(  RUNNING_MEDIANS[R,COL_L] - data[R,COL_L]  )
		}

	################################
	## MEDIAN ABSOLUTE DEVIATIONS ##
	################################
	MEDIAN_ABSOLUTE_DEVIATIONS = matrix( NA, ncol=NCOL, nrow=NROW, dimnames=dimnames(data))
 	for( L in LAGS ) 
		for( R in (HISTORY_DAYS+1):NROW ) {
			COL_L = paste(L, sep="")
			TEMP = ABSOLUTE_DEVIATIONS[(R-HISTORY_DAYS):(R-1), COL_L]
			if(  sum(!is.na( TEMP ) ) >= HISTORY_DAYS/2  ) 
				MEDIAN_ABSOLUTE_DEVIATIONS[R, COL_L] = median( TEMP, na.rm=T)
		}

	##################################
	## Plotting the running medians ##
	##################################
	Y_MIN     = min( as.vector(data), na.rm=T)
	Y_MAX     = max( as.vector(data), na.rm=T)
	Y_MAX_MAD = max( as.vector(MEDIAN_ABSOLUTE_DEVIATIONS), na.rm=T)
	Y_MAX 	  = Y_MAX + 2*Y_MAX_MAD
	X 		  = 1:NROW
	ZEROES = rep(0, times=length(X))		 

	## x-axis labels.
	X_NUMBER_OF_LABELS = 10
	X_JUMP = round( (NROW/X_NUMBER_OF_LABELS)/10, 0 )*10
	X_TICK_PLACES	 = (0:(X_NUMBER_OF_LABELS)) * X_JUMP
	X_LABELS = X_TICK_PLACES	
	X_LABEL_ORIENTATION = 1


	## y-axis labels.
	CEILING = ceiling(Y_MAX/1000)*1000
	Y_JUMP = round(CEILING/5, 0)
	if( Y_JUMP == 0 ) Y_JUMP = 500
	## These are the plot settings if "horizontal==T" and there are more than 10 lags.
	Y_TICK_PLACES = seq(0, CEILING, by=Y_JUMP)
	Y_LABELS = Y_TICK_PLACES
	Y_FONT_AXIS = 1

	## Reset x-axis labels if start date is provided.
	START_DATE = accrued_data[["start"]][[1]]
	if( START_DATE != 1 ){
			# x-axis labels need to be dates.
			X_LABELS 	= as.Date(START_DATE + X_TICK_PLACES, origin=as.Date("1970-01-01")) 
			X_LABEL_ORIENTATION = 2
	}		


	for( L in LAGS ) {
		COL_L = paste(L, sep="")

		if( accrued_data[["start"]] != 1 ) par(mar=c(5.0,2.5,2,1) )
		else par(mar=c(2.0,2.5,2,1) )
		

		## This asks user if (s)he is ready to view the next plot.
		if( names(dev.cur()) %in% deviceIsInteractive() ) devAskNewPage(ask = TRUE)

		plot( X, data[,COL_L],  
		      xlab="", ylab="Value", 
		      ylim=c(Y_MIN,Y_MAX), main=paste("Lag=",L,sep=""), axes=FALSE, 
			  frame.plot=TRUE, cex.main=1, pch=16, cex=0.5, col="dimgrey", ...)		 
		lines(X, RUNNING_MEDIANS[,COL_L], col = "blue", lwd=2)
		lines(X, pmax(ZEROES,RUNNING_MEDIANS[,COL_L]-2*MEDIAN_ABSOLUTE_DEVIATIONS[,COL_L]), col="darkgreen", lty=1)
		lines(X, RUNNING_MEDIANS[,COL_L]+2*MEDIAN_ABSOLUTE_DEVIATIONS[,COL_L], col="darkgreen", lty=1)
		abline(a=0, b=0)
		axis( 1, at=X_TICK_PLACES, labels=X_LABELS, las=X_LABEL_ORIENTATION, font.axis=1, cex=0.8)
		axis( 2, at=Y_TICK_PLACES, labels=Y_LABELS, font.axis=Y_FONT_AXIS, las=2)

	} # END for( L in LAGS )
	
	if( names(dev.cur()) %in% deviceIsInteractive() ) devAskNewPage(ask = FALSE)

}

