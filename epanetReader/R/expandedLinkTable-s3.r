#************************************
#
#  (C) Copyright IBM Corp. 2015
#
#  Author: Bradley J Eck
#
#************************************

#  File:  expandedLinkTable-s3.r
#
#  Purpose: define an s3 class for a link table 
#           with coordinates added


#' Expanded Link Table
#' 
#' Create an expandedLinkTable object by adding node coordinates to a 
#' data frame of pipes, pumps, or valves.  
#'
#' @export
#' @param Links data frame of Pipes, Pumps or Valves of from epanet.inp 
#' @param Coordinates table of epanet.inp
#' @return an expandedLinkTable object 
#' @examples 
#' x <- expandedLinkTable(Net1$Pipes, Net1$Coordinates) 
#' print(x)
#' plot(x) 
expandedLinkTable <- function( Links, Coordinates ){

  # handle a missing table 
  if( is.null(Links) ) {
    return(NA)
	#ept <- NA
  } else {
    ept <- merge( x = Links, by.x = "Node1", all.x = TRUE, sort = FALSE, 
                  y = Coordinates, by.y = "Node" ) 
    #rename
    names(ept)[ grep("X.coord", names(ept)) ]  <- "x1"
    names(ept)[ grep("Y.coord", names(ept)) ]  <- "y1"
    
    #Node2 coords
    ept <- merge( x = ept, by.x = "Node2", all.x = TRUE, sort = FALSE, 
                  y = Coordinates, by.y = "Node" ) 
    #rename
    names(ept)[ grep("X.coord", names(ept)) ]  <- "x2"
    names(ept)[ grep("Y.coord", names(ept)) ]  <- "y2"
    
    # midpoints for labeling 
    ept$midx <- (ept$x1 + ept$x2) / 2 
    ept$midy <- (ept$y1 + ept$y2) / 2  
	
	# put the columns into order 
	ept <- ept[ ,c( names(Links), 'x1', 'y1', 'x2', 'y2', 'midx', 'midy') ]
    
  }
      
    class(ept) <- c("expandedLinkTable", "data.frame")
    
    return(ept)
}

#' plot an expanded link table 
#' 
#' @export
#' @param x object of type expandedLinkTable
#' @param add logical indicating whether to add to the currently active plot.  
#'        add=FALSE creates a new plot.
#' @param label logical indicating if the links should be labeled at the mid points
#' @param linewidths passed to lwd argument in segments()
#' @param ... further arguments passed to segments() 
#' @details 
#' An implementation of the generic plot function for
#' expandedLinkTable objects. Links are drawn using segments(). 
#' Useful for building up network plots.
plot.expandedLinkTable <- function(x, add=FALSE, label=FALSE, linewidths = 3, ...){
  
    if( add == FALSE ){
      # generate a blank plot first 
      graphics::plot( range( c(x$x1, x$x2) ),
            range( c(x$y1, x$y2) ),
            type = 'n',
            xlab = "", xaxt = 'n',
            ylab = "", yaxt = 'n'
			)
      
    } 
    
    # just put the segments out there 
    graphics::segments( x0 = x$x1, y0 = x$y1,
              x1 = x$x2, y1 = x$y2,
			  lwd = linewidths, ... )  
               
	
    if( label == TRUE ){
      graphics::text( x$midx, x$midy, x$ID)
    }
}
