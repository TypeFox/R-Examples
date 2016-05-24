#************************************
#
#  (C) Copyright IBM Corp. 2015
#
#  Author: Bradley J Eck
#
#************************************

# constructor for s3 class  
#' Sparkline Table 
#' 
#' Create S3 object of data for table of sparklines 
#' 
#' @export 
#' @param df data.frame of values to plot.   
#' @param row.var variable for rows of the table 
#' @param col.vars variables for columns of the table 
#' @param xvar optional name of variable for horizontal axis of sparkline plots
#' @param xrange.labels optional vector of length 2 with labels for the first
#'        and last quantities plotted on x-axis, often a date and/or time
#' @seealso plotSparklineTable
#' 
sparklineTable <- function( df, row.var, col.vars, xvar = NULL, xrange.labels = NULL ){
	
	#argument checking
	if( is.data.frame(df) == FALSE) stop("df must be a data.frame")
	if( length( row.var) != 1 ) stop("row.var must have length 1")
	if( (row.var %in% names(df))  == FALSE ) stop(paste("row.var", row.var, "is not a column in data.frame df"))
	if( min( col.vars %in% names(df)) < 1 ) stop(paste("at least one value in col.vars is not a column in data.frame df")) 
	if( is.null(xvar) == FALSE ){ 
		if( (xvar %in% names(df))  == FALSE ) {
			stop(paste("xvar", xvar, "is not a column in data.frame df"))
		}
	}
	# x range label  
	xrl <- c("","")
	if( is.null(xrange.labels) == FALSE ){
		if( length(xrange.labels) != 2 ) stop("xrange.labels must have length 2")
		xrl <- as.character(xrange.labels)
	}
	
	
	urv <- unique( df[ , row.var])
	Nrv <- length(urv ) 	
	Ncv <- length(col.vars)
	M <- getLayoutMatrix( num.row.var = Nrv,
			num.col.vars = Ncv  )
	# number of sparklines 
	Nsl <- Nrv * Ncv
	
	# need to create data structure for
	# sparkline table
	
	###################
	# Constructor 
	###################
	
	# list of sparkline matrices 
	sll <- list()
	
	k = 1 
	# Loop over the rows 
	for( i in 1:Nrv){
		# Loop thru the params 
		for(j in 1:Ncv){
			# first param, start value 
			xy <- sparkline( df, row.var, urv[i], col.vars[j], xvar)
			
			sll[[k]] <- xy
			k <- k + 1 
		}
	}
	
	# check the data 
	sparklineDataCheck( sll )

	# make the table header 
	H <-  row.var
	for( j in 1:Ncv){
	
	    H <- c(H, xrl[1])
		H <- c(H, col.vars[j])
	    H <- c(H, xrl[2])	
		
		#plotWord( xrl[1], font=2)
		#plotWord( col.vars[j], font = 2 , cex = 1)
		#plotWord( xrl[2], font=2)
	}
	
    # build up the object to return 	
	slt <- list( sparklines = sll,  
			     layoutMatrix = M,
				 tableHeader = H, 
				 rowLabels = urv)
	
	class( slt ) <- 'sparklineTable'

	return( slt)
}

#' Plot Sparkline Table 
#' 
#' @export 
#' @param x object of class sparklineTable 
#' @param ... further arguments passed to par
plot.sparklineTable <- function( x,... ){

	#########################
	#   Plotting 
	#########################
	# create the plot grid 
	oldpar <- graphics::par(no.readonly = TRUE) # keep up w the old params 
	graphics::par( mar = c(0,0,0,0), oma = rep(1,4),...)
	graphics::layout( mat = x$layoutMatrix, respect = FALSE  )
	
	lapply( x$tableHeader, plotWord, font = 2)
	
	
	# Loop over the rows 
	Nrv <- length( x$rowLabels)
	Ncv <- length(x$sparklines) / Nrv 
    k = 1 	
	for( i in 1:Nrv){
		
		plotWord( x$rowLabels[i]) 
		
		# Loop thru the params 
		for(j in 1:Ncv){
			# first param, start value 
			xy <- x$sparklines[[k]] 
			N <- dim(xy)[1]
			plotWord( xy[1,2] ) 
			graphics::plot( xy )
			plotWord( xy[N,2] ) 	
			k = k + 1 
		}
	}
	
	# go back to the old params 
    graphics::par(oldpar)
	
}

	
	
sparklineDataCheck <- function( slt ){
# check data for which sparklines will be plotted to make 
# sure visualization will be faithful 
	
	#  (1) Do all sparklines have the same range of xvar ? 
	kmax <- length( slt )
	
	# calc the x range for each sparkline   
	xrange <- list()	
	for( k in 1:kmax){
		xrange[[k]] <- range( slt[[k]][,1] )
	}			
	
	# see if they are all the same 
	sameXrange <- ( length( unique(xrange)) == 1 ) 
	
	if( sameXrange == FALSE ){
		warning( "Sparklines have different ranges of X") 
	}
}
	
getLayoutMatrix <- function( num.row.var, num.col.vars ){
	
	
	# determine the matrix size for the layout 
	nr <-  1 + num.row.var  
	nc <- 1 + 3 * num.col.vars
	
	# col labels spread over three
	#header <- c(1, unlist( lapply ( 2:(num.col.vars+1), FUN = rep, times = 3 ) ) ) # param name takes up 3 columns 
	#ki2j1 <- tail( header, 1 ) + 1  
	#vals <- c(header, seq( from = ki2j1, by = 1, length.out = (nr-1) * nc ) )  
	
	vals <- 1 : (nr * nc)	
	
	M <- matrix( data =vals, nrow = nr, ncol = nc, byrow = TRUE) 
	
	return( M )
}

plotWord <- function(w, ...){
	
	if( is.numeric(w) ) {
		
		w <- format( w, digits = 3)
	} 
	
    graphics::plot(c(0,1),c(0,1), type = 'n', 
			xaxt='n', yaxt='n', xlab = '', ylab = '', 
			frame.plot = FALSE )
	graphics::text( .5, .5, w, ...) 
	
}	
	
