#************************************
#
#  (C) Copyright IBM Corp. 2015
#
#  Author: Bradley J Eck
#
#************************************

#
#  File:  epanet.inp-s3.r
#
#  Purpose:  S3 implemetation of epanet.inp



#' Read .inp file 
#' 
#' Read an Epanet .inp file into R
#'
#' @aliases epanet.inp
#' @export
#' @param file the name of the file to read 
#' @details 
#' This function reads a text file in Epanet's .inp format 
#' and returns an S3 object with entries for 
#' sections of the .inp file.  Sections of the .inp file that are implemented 
#' appear in the Value section.
#' 
#' Fields for node or link ID are stored as characters not factors or integers.
#' However, some fields are stored as factors
#' to allow more informative summaries. Examples include valve type and pipe status.  
#'        
#' Sections that are absent from the .inp file are NULL in the list.
#' 
#' Columns of data.frames use the headings exported by 
#' the Epanet GUI. 
#' 
#' The [OPTIONS] section in the .inp file is used to update
#'          a list of Epanet's default options. In this way if an option such as
#'          units is not specified by the .inp file, the units that would be used by
#'          default are provided.
#' 
#' In the [PATTERNS] and [CURVES] sections, integers used as names of list elements are backquoted
#'          according to the default behavior in R.  So if the .inp file has a pattern "1"
#'          this pattern will appear as element `1` in the list that is returned. A warning is issued in this case. 
#' 
#' @return Returns an epanet.inp S3 object with 
#' elements of the following names and types corresponding to sections
#' of the .inp file. Sections missing from the .inp file have a value of NULL.
#' \item{Title}{character}
#' \item{Junctions}{data.frame}
#' \item{Tanks}{data.frame}
#' \item{Reservoirs}{data.frame}
#' \item{Pipes}{data.frame}
#' \item{Pumps}{data.frame}
#' \item{Valves}{data.frame}
#' \item{Demands}{data.frame}
#' \item{Status}{data.frame}
#' \item{Patterns}{list}
#' \item{Curves}{list}
#' \item{Energy}{character}
#' \item{Times}{character}
#' \item{Options}{list}
#' \item{Coordinates}{data.frame}
#' 
#' @references Rossman, L. A. (2000). Epanet 2 users manual. US EPA, Cincinnati, Ohio.
#' 
#' http://nepis.epa.gov/Adobe/PDF/P1007WWU.pdf
#' 
#' @examples
#' # path to Net1.inp example file included with this package
#' inp <- file.path( find.package("epanetReader"), "extdata","Net1.inp") 
#' 
#' #read the network file into R
#' n1 <- read.inp(inp)
#' summary(n1)
#' names(n1)
#' summary(n1$Junctions)
#' summary(n1$Pipes)
#' plot(n1) 
read.inp <- function (file ){
 
  return(  epanet.inp( file ) )
  
}

epanet.inp <- function( file ){
	
	allLines <- readLines( file )
  
  # read in all the sections  
  titl <- TITLE( allLines )
  junc <- JUNCTIONS( allLines )
  tank <- TANKS (allLines)
  resr <- RESERVOIRS( allLines) 
  pipe <- PIPES(allLines)
  pump <- PUMPS(allLines)
  valv <- VALVES(allLines)
  dmd  <- DEMANDS(allLines) 
  pats <- PATTERNS(allLines)
  crvs <- CURVES(allLines)
  #ctrl
  engy <- ENERGY(allLines)
  stat <- STATUS(allLines)
  #emit
  #qlty
  #srcs
  #rxns
  #mix
  tims <- TIMES(allLines)
  #rpt
  opts <- OPTIONS(allLines)
  coor <- COORDINATES(allLines)
  #vert
  #labs
  #bdrp
  
  # make a list of all the elements 
  inp <- list( Title = titl,
               Junctions = junc,
               Tanks = tank,
               Reservoirs = resr,
               Pipes = pipe,
               Pumps = pump,
               Valves = valv,
			   Demands = dmd, 
               Patterns = pats,
               Curves = crvs,
               Energy = engy,
			   Status = stat,
               Times = tims,
               Options = opts,
               Coordinates = coor)               
  
  class(inp) <- "epanet.inp"
  return( inp )
}


#' Summary Method for epanet.inp
#'
#' Summarizes the network by printing the Title of the
#' network and the number of each type of elements.
#'
#' @export
#' @param object of class epanet.inp
#' @param ... futher arguments
summary.epanet.inp <- function( object, ... ){
 
  Number <- unlist(lapply( object, function(object) dim(object)[1]))
  counts <- as.data.frame(Number)
  
   inpSmry <- list( Title = object$Title,
				   entryCounts = counts )
	class(inpSmry) <- "summary.epanet.inp"	   
	
	return(inpSmry)
	
}

print.summary.epanet.inp <- function(object, ...){
	
  cat( object$Title, sep = "\n")
  cat( "\n")
  
  print( object$entryCounts)
}

#' Plot .inp Links
#'
#' Add lines for pipes, pumps and valves
#'  from an epanet.inp object to an existing plot 
#'  
#' @export 
#' @param x epanet.inp object 
#' @details 
#' Helper function for building up a plot of the network by
#' adding links to an existing plot.  
#' @examples 
#' ## make a new blank plot 
#' plot( range(Net1$Coordinates$X), range(Net1$Coordinates$Y), type = 'n') 
#' ## add the links
#' plotInpLinks(Net1) 
plotInpLinks <- function(x){
	
  #############
  #  Pipes  
  ############# 
  if( is.null( x$Pipes) == FALSE ){
    graphics::plot( expandedLinkTable( x$Pipes, x$Coordinates ),
          add = TRUE,  
          label = FALSE ) 
  }

  #############
  #  Pumps  
  ############# 
  if( is.null( x$Pumps )  == FALSE ){
	  ept <-  expandedLinkTable( x$Pumps, x$Coordinates )
	  graphics::plot( ept,
			  add = TRUE,  
			  label = FALSE) 
	  graphics::points( ept$midx, ept$midy, pch = 8 ) 
  }

  
  #############
  #  Valves 
  ############# 
  if( is.null( x$Valves )   == FALSE ){
	  evt <- expandedLinkTable(x$Valves, x$Coordinates)
	  graphics::plot( evt,
			  add = TRUE,  
			  label = FALSE) 
	  
	  graphics::points( evt$midx, evt$midy, pch = 25 ,
			  bg="black", col = "black" )  
  }
}

#' Plot Node Elements 
#'
#' Adds node elements from epanet.inp object to an existing plot
#' 
#' @export 
#' @param x epanet.inp object
#' @param plot.junctions logical indicating whether to plot junctions 
#' @details  Helper function for building up a network plot. Tanks and
#' Reservoirs are shown using plot characters (pch) '           16 and 15.
#' Junctions, if plotted, appear as pch ="."
#' @examples 
#' ## make a new blank plot 
#' plot( range(Net1$Coordinates$X), range(Net1$Coordinates$Y), type = 'n') 
#' ## add the nodes, including junctions 
#' plotInpNodes(Net1, TRUE ) 
plotInpNodes <- function( x, plot.junctions){
  ############# 
  #  Junctions 
  #############
  if( plot.junctions ){
     
    jpts <- merge( x = x$Junctions, by.x = "ID", all.x = TRUE,
                   y = x$Coordinates, by.y = "Node" )

    graphics::points( jpts$X.coord, jpts$Y.coord, pch = "." )

  }

  ###########
  #  Tanks 
  ###########
  if( is.null ( x$Tanks ) == FALSE ){
    tpts <- merge( x = x$Tanks, by.x = "ID", all.x = TRUE,
                   y = x$Coordinates, by.y = "Node" )
    graphics::points( tpts$X.coord, tpts$Y.coord, pch = 16, col='black' )
  }

  ##############
  #  Reservoirs
  ##############
  if( is.null( x$Reservoirs) == FALSE ){
    rpts <- merge( x = x$Reservoirs, by.x = "ID", all.x = TRUE,
                   y = x$Coordinates, by.y = "Node" )
    graphics::points( rpts$X.coord, rpts$Y.coord, pch = 15, col='black' )
  }
}


#' Plot Legend of Network Elements
#' 
#' Add legend of network elements to the active plot 
#' 
#' @export  
#' @param legend.locn keyword for location of legend. See details of legend()
#'        function.
#' @details
#' Helper function for adding a legend to the active plot.  
#' Uses plot characters 16, 15, 8 and 25 for Tanks, Reservoirs, Pumps and Valves 
#' for compatibility with plotInpNodes() 
#' @examples
#' ## make a new blank plot 
#' plot( c(0,1), c(0,1), type = 'n') 
#' ## add the nodes, including junctions 
#' plotElementsLegend('topright') 
plotElementsLegend <- function(legend.locn) {
	
  graphics::legend( legend.locn, bty = 'n',
          c("Tanks", "Reservoirs", "Pumps", "Valves"),
         pch = c(16, 15, 8, 25),
         pt.bg = 'black',
         col = 'black')
}

#' Plot Method for epanet.inp 
#' 
#' Make a plot of the network using base graphics
#' 
#' @export
#' @param x object of class epanet.inp 
#' @param plot.junctions logical indicating whether to plot junctions
#' @param legend.locn character string passed to legend() specifying
#'        the location of the legend on the plot 
#' @param ... other arguments passed to plot()
#' @details
#' Implements the generic plot function for S3 objects of class epanet.inp.
#' The plot is built from base graphics by creating a blank plot and then calling 
#' the helper functions plotInpLinks(), plotInpNodes(), plotElementsLegend().  
#' @examples
#' plot(Net1) 
plot.epanet.inp <- function( x, 
                             plot.junctions  = TRUE,
                             legend.locn = "topright",
                                ... ) {

  
  # check that coordinates actually exist   
  if( is.null(dim(x$Coordinates) )){
    stop("network does not have coordinates" )
  }
  
  ### Do we have vertices ? 
  

  
  # create blank plot 
  graphics::plot( range(x$Coordinates$X),
        range(x$Coordinates$Y),
        type = 'n', asp = 1, 
        xlab = "", xaxt = 'n',
        ylab = "", yaxt = 'n', ... )

  plotInpLinks(x)
  
  plotInpNodes(x, plot.junctions)

  plotElementsLegend(legend.locn) 
  
}



