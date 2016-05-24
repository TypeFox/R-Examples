#************************************
#
#  (C) Copyright IBM Corp. 2015
#
#  Author: Bradley J Eck
#
#************************************

#'  Sparkline
#' 
#' Create sparkline object by extracting from a data frame
#' 
#' @export 
#' @param df data.frame from which data for the sparkline is extracted    
#' @param id.var variable in df with IDs 
#' @param ID value in id.var on which to extract 
#' @param yvar name of variable for the y values in the sparkline  
#' @param xvar optional name of variable for horizontal axis of sparkline plots
#' @details 
#' Creates an object with info for a single sparkline by extracting 
#' from a data.frame.  The function works on data.frames with one column of ID variables
#'          and possibly several columns of other variables.  The main use is 
#'          as a helper function for building up a \link{sparklineTable}. 
#' 
#' @examples 
#' ## look at the names in the built-in data set Theoph
#' names(Theoph) 
#' ## make sparkline object for the concentration over time in subject 2
#' sl <- sparkline(df= Theoph, id.var = 'Subject', ID = 2, yvar='conc', xvar = 'Time') 
#' plot(sl)
sparkline <- function( df, id.var, ID, yvar, xvar){
	
	D <- df[ which( df[ , id.var] == ID ) , c( xvar, yvar )   ] 
	D <- as.matrix(D)	
	# if D has only one column we assume it's already in the right order
	# and add a vector of indices to plot against 
	# if D has two columns we sort it into xvar order
	
	Dncol <- dim(D)[2]
	Dnrow <- dim(D)[1]
	
	colvarname <- paste(ID, yvar)
	
	if( Dncol == 1 ){ 
		
		D <- cbind( 1:Dnrow, D)
		
		colnames(D) <- c("index", colvarname) 
		
	} else if( Dncol == 2) { 
		
		xord <- order( D[ , 1])
		D <- D[xord, ]
		
		colnames(D) <- c(xvar, colvarname)
		
	} else { 
		stop("D should only have two cols, something is wrong")
	}
	
	class(D) <- 'sparkline'
	
	return(D)
	
}

#' Plot a sparkline 
#' 
#' @export 
#' @param x sparkline object   
#' @param ... further arguments passed to plot.default 
#' @details 
#' Implementation of the generic plot function for a single sparkline object.
#' The primarily used to build up plots of a sparklineTable
#' @seealso sparkline
plot.sparkline <- function( x, ... ){	
	
	dimxy <- dim(x)
	N <- dimxy[1]
    # argument checking 	
	if( dimxy[2] != 2 ) stop(" input x must have two columns")
	
	# give a 10% buffer in the vertical 
	yr <- range(x[,2])
	ybuff <- ( yr[2] - yr[1] ) * 0.10 
	ylimits <- c( yr[1] - ybuff, yr[2] + ybuff )
	
	# give a 1% buffer in the horizontal 
	xr <- range(x[,1])
	xbuff <- (xr[2] - xr[1]) * 0.01
	xlimits <- c(xr[1] - xbuff, xr[2] + xbuff ) 
	
	graphics::plot.default(x, type = 'l', col = 'gray',  
			xaxt='n', yaxt='n', xlab = '', ylab = '', 
			ylim =  ylimits, xlim = xlimits,
			frame.plot = FALSE, ... )
	graphics::points( x[1,1], x[1,2], pch = 16, cex = .9 )
	graphics::points( x[N,1], x[N,2], pch = 16, cex = .9 )  
}

