lagHistogram <- function(x, maxLagPlot=NULL, ...) {

	## Throw an error if the argument is not of the correct class.
	if( class(x) != "accrued" )  stop("ERROR: argument is not an object of the 'accrued' class.")

	accrued_data = x
	
	data = accrued_data[["data"]]
	final = accrued_data[["final"]] 
	MAX_LAGS = ncol(data) - 1
	NCOL = ncol(data)
	NROW = nrow(data)

	## Each histogram has 10 bins: [0,10), [10,20), [20,30), [30,40)...,[90,100].
	## We populate the bins for every lag, however, not all lags will be plotted.
	NBINS = 10 
	HistogramBinData = matrix( data=0, nrow=(NBINS+1), ncol=NCOL, dimnames=list(1:(NBINS+1),0:MAX_LAGS) )


	for( L in 1:(MAX_LAGS+1) ) {
	
		laggedData = data[,L]

		for( R in 1:NROW ) {
	
			PROPORTION = NA
	
			## If there are final counts or proxy final counts for that date.	
			if(  !is.na( final[R] ) & ( final[R] > 0 )  ) { 
	
				## The last row counts the number of non-NA data points in the histogram.
				HistogramBinData[ NBINS+1, L ] = HistogramBinData[ NBINS+1, L ] + 1
				PROPORTION = laggedData[R] / final[R]
	
				## If there is data for that encounter date, but "NA" for that lag, treat the NA as 0.				
				if( is.na(PROPORTION) )  PROPORTION = 0

				## Populate the HistogramBinData structure, using the [,) convention for histograms.
				for( B in 1:NBINS ) {
					if( B == 1 ){
						if( PROPORTION < B/NBINS )  HistogramBinData[B, L] = HistogramBinData[B, L] + 1		
					} else if ( B < 10 ) {
						if( ( PROPORTION >= (B-1)/NBINS ) & ( PROPORTION < B/NBINS )  )  
							HistogramBinData[B, L] = HistogramBinData[B, L] + 1
					} else {
					# Case: B = 10
						if( ( PROPORTION >= (B-1)/NBINS ) & ( PROPORTION <= B/NBINS )  )  
					 		HistogramBinData[B, L] = HistogramBinData[B, L] + 1
					} # END  if-else (B==1) 
				} # END  for( B in 1:NBINS )

			} # END  if(  !is.na( final[R] ) & ( final[R] > 0 )  ) 
		} # END  for( R in 1:NROW ) )
	} # END  	for( L in 0:MAX_LAGS )


	
	## Data structure has been created.
	## Now set the graphical parameters.

	NDAYS = min(c(14, MAX_LAGS))
	if( !is.null(maxLagPlot) )  NDAYS = max(1, maxLagPlot)
	NROW = dim(HistogramBinData)[[1]]
	NBINS = NROW-1
		
	#############################################################################################################
	# Grid settings.
	#############################################################################################################
	# Width of left column (used for lags).
	LEFT_WIDTH = 0.1

	# Width of right column (used for histograms and axes).
	# A little bit of right-hand padding is there so the right-most
	# axis labels doesn't get cut off. 
	RIGHT_WIDTH = 1 - LEFT_WIDTH*(1.2)

	# Height of bottom rectangle  
	BOTTOM_HEIGHT = 0.15
	if( NDAYS < 14 )  BOTTOM_HEIGHT = 0.15+(14-NDAYS)/120


	# Bottom Axis Margin shift--this is hard to optimize. 
	AXIS_MARGIN = 2.25


	# Height of upper rectangle
	TOP_HEIGHT = 1-BOTTOM_HEIGHT
	
	# Increment from histogram to histogram
	HIST_JUMP   = 1/( NDAYS + 0.5 )
	if( NDAYS < 4 ) HIST_JUMP = min( HIST_JUMP, 1/3 )


	# Height of individual histograms --- This may change.
	# SHRINK CAN BE SET TO OTHER VALUES, SUCH AS 0.95
	SHRINK = 1
	HIST_HEIGHT = SHRINK*HIST_JUMP 
	
	### Clear the current plotting region so grid does not just overlay if a plot already exists.
	### Is there not a better way to do this?
	plot.new() 
	plot.window(c(0,1),c(-1,1),mar = c(1,0,0,0) )

	#############################################################################################################
	### Start grid.
	#############################################################################################################
	grid.rect( gp=gpar(col="white") )
	
	## NOTE THAT BY CHOOSING A JUSTIFICATION OF "left" "bottom"
	## GIVES COORDINATES WITHIN EACH VIEWPORT THAT COINCIDE WITH THE 1ST 
	## QUADRANT OF THE CARTESIAN PLANE, AND ARE ALWAYS RELATIVE TO THE
	## VIEWPORT THAT IS THE NEXT LEVEL UP.

	#############################################################################################################
	# Outer Rectangle (with black border on the outside)
	#############################################################################################################
	## This is for aesthetics
	pushViewport(  viewport( x=0.01, y=0.01, width=0.98, height=0.98, just=c("left","bottom") )  )
		grid.rect( gp=gpar(col="black") )

		#########################################################################################################
		# Inner Rectangle (no border, buffer between outer border and inner border)
		#########################################################################################################
		pushViewport(  viewport( x=0.01, y=0.01, width=0.98, height=0.98, just=c("left","bottom") )  )


			#####################################################################################################
			# Lower Portion of Inner Rectangle (SITE ID and X-AXIS -- BIN SIZES)
			# Coordinates are elative to the "inner rectanlge", since it's within IR's viewport.
			#####################################################################################################
			pushViewport(  viewport( x=0, y=0, width=1, height=BOTTOM_HEIGHT, just=c("left","bottom") )  )

				#################################################################################################
				## Right Column of Lower Rectangle: X-AXIS & LABELS 
				#################################################################################################
				pushViewport(  viewport( x=LEFT_WIDTH, y=0, width=RIGHT_WIDTH, height=1, just=c("left","bottom") )  )
					pushViewport( plotViewport( c(AXIS_MARGIN, 0, 0, 0) )   )  
						FONT_SIZE = 10
						#if( NDAYS < 14 ) { FONT_SIZE = 5 }
						grid.xaxis( label=T, at=seq(0, 1, by=1/NBINS), gp=gpar(fontsize=10) )
					upViewport() ## End Bin Axis
				upViewport() ## End of Right Column of Lower Rectangle
			upViewport() ## End of Lower Rectangle

	
			#####################################################################################################
			# Upper Portion of inner rectangle (histograms)
			#####################################################################################################
			## Coordinates are relative to the inner rectangle, and since the lower rectangle 
			pushViewport(  viewport( x=0, y=BOTTOM_HEIGHT, width=1, height=TOP_HEIGHT, just=c("left","bottom") )  )
	
				#################################################################################################
				## Histogram Rows of Upper Portion.
				#################################################################################################
				for( L in 0:(NDAYS+1) ) {
					if( L <= NDAYS ) {  
						NAME = as.character(L)
						#########################################################################################
						## Left-hand column of Histogram Rows (for lag labels), surrounded by boxes.
						#########################################################################################
						pushViewport(  viewport( x=0, y=(L-0.75)*HIST_JUMP, width=LEFT_WIDTH, height=HIST_HEIGHT, 
									  just=c("left","bottom") )  )
							grid.rect(  gp=gpar( col="forestgreen")  )
							grid.text( NAME )
						upViewport()
					} # END  if( L <= NDAYS )

					##############################################################################################
					## Right-hand of Histogram Rows (for the histograms)
					##############################################################################################
					pushViewport(  viewport( x=LEFT_WIDTH, y=(L-0.75)*HIST_JUMP, width=RIGHT_WIDTH, height=HIST_HEIGHT, 
								  just=c("left","bottom") )  )
						for( B in 1:NBINS ) {
							######################################################################################
							## Histogram Cell for lag L and bin B
							######################################################################################
							pushViewport(  viewport( x=(B-1)/NBINS, y=0, width=1/NBINS, height=1, just=c("left","bottom") )  )
								if( L <= NDAYS ) { 
									##############################################################################
									## Inner Viewport where Cell is Drawn
									## Height is given by the proportion.
									##############################################################################

									BAR_HEIGHT = ( HistogramBinData[B, L+1] / HistogramBinData[NROW, L+1] )
									if( is.na(BAR_HEIGHT) ) BAR_HEIGHT = 0
									pushViewport(  viewport( x=0, y=0, width=1, height=BAR_HEIGHT, just=c("left","bottom") )  )
										COLOR = NULL
										if( HistogramBinData[B, L+1] > 0 ) COLOR = "dodgerblue"
										else COLOR = "grey"	
										grid.rect( gp=gpar(fill=COLOR, col=COLOR) )
									upViewport()

									###############################################################################
									## Viewport for lower gray line, so zero of highest bin is visible, even when
									## data are present.
									###############################################################################
									pushViewport(  viewport( x=0, y=0, width=1, height=0.01, just=c("left","bottom") )  )
										COLOR = "grey"
										grid.rect(  gp=gpar( fill=COLOR, col=COLOR )  )
									upViewport()										

								} else {
									pushViewport(  viewport( x=0, y=0, width=1, height=0.01, just=c("left","bottom") )  )
										COLOR = "grey" 
										grid.rect( gp=gpar(fill=COLOR, col=COLOR) )
									upViewport()
								} # END  if-else( L <= NDAYS ) 			

							upViewport()
							
						} # END  for( B in 1:NBINS ) 
						
					upViewport()
				} # END  for( L in 0:(NDAYS+1) )

			upViewport() ## End of upper rectangle

		upViewport() ## End of inner rectangle
	upViewport() ## End of outer rectangle

	#dev.off()  ## Only used if the pdf part is used.


}
