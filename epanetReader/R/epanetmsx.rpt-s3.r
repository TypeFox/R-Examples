#************************************
#
#  (C) Copyright IBM Corp. 2015
#
#  Author: Bradley J Eck
#
#************************************

#' Read msx results 
#' 
#' reads an Epanet-msx .rpt file into R
#' 
#' @aliases epanetmsx.rpt 
#' @export 
#' @param file the name of the file to read
#' @return Returns an epanetmsx.rpt S3 object . 
#'
#' \item{nodeResults}{data.frame}
#' \item{linkResults}{data.frame}
#' @details  Specify the needed outputs from an Epanet-msx simulation in the
#'  [REPORT] section of the .msx file to create reports for reading with 
#' with this function. 
#'
#' The function returns an S3 object (list) with a data.frame for node results and
#' data.frame for link results.  These data.frames contain results from all
#' the time periods to facilitate time series plots. 
#'
#' @references Shang, F., Uber, J.G., Rossman, L.A. (2011)
#'             EPANET Multi-species Extension User's Manual.
#'             US Environmental Protection Agency, Cincinnati. 
#' @examples
#' # path to example file included with this package
#' msr <- file.path( find.package("epanetReader"), "extdata","example.rpt") 
#' 
#' #read the results into R
#' x <- read.msxrpt(msr)
#' names(x)
#' summary(x)
#' plot(x)
read.msxrpt <-function( file ){
	
	mro <- epanetmsx.rpt(file)
	return( mro)
	
}


epanetmsx.rpt <- function( file ) {
# this based heavily on the epanet.rpt function
	
  # read all the lines in the file 
  allLines <- readLines(file)
  #checkRptFile( allLines ) 

  title <- getTitle( allLines) 

  lengthOfAllLines <- length( allLines)
  
  resLines <- grep("<<<.+>>>", allLines)
  numTables <- length(resLines)

  # empty lists of results 
  nodeResList <- list()  
  linkResList <- list()  
  
    # initialize indices for these lists 
  ni <- 0  # node index
  li <- 0  # link index
  
  # go through the tables  
  for( i in 1:numTables ){
    # get the section  
    sectRange <- .getSectionRange( i, resLines, lengthOfAllLines)
    sect <- allLines[ sectRange$start : sectRange$end ]  
    
    # create a data frame from this section
    df <- msxSection2df( sect )  

    # decide if it's for link or node results 
    isNODE <- grepl( "Node", sect[1] )
    isLINK <- grepl( "Link", sect[1] )
    
    #add  data to approriate list of data frames 
    if( isNODE ){
      # increment node indexer 
      ni <- ni + 1
      nodeResList[[ ni ]] <- df
    }
    
    if( isLINK ){
      # increment indexer 
      li <- li + 1
      linkResList[[ li ]] <- df
    }
  }
  
  # combine all these results together 
  nodeResDF <- do.call("rbind", nodeResList )
  linkResDF <- do.call("rbind", linkResList )
 
  # and make a list of the to return 
  allResults <- list( Title = title, 
                      nodeResults = nodeResDF, 
                      linkResults = linkResDF )
  
  class(allResults) <- "epanetmsx.rpt"

  return( allResults ) 
}


#' Summary of Epanet-msx Simulation Results
#'
#' Provides a basic summary of simulation results 
#'
#' @export
#' @param  object of epanetmsx.rpt class
#' @param ... further arguments passed to summary()
summary.epanetmsx.rpt <- function( object, ...){

  # time info  (you have results from tmin to tmax at an interval of delta_t)



  ################
  #  node results
  ################
  
  # unique node IDs 
  uni <- unique(object$nodeResults$ID)
  if( length(uni) > 0 ){
	  # node result names
	  nrn <- names(object$nodeResults) 
	  
	  # node result summary over species  
	  jmax <- length(nrn) - 1 
	  # species results
	  nsr <- object$nodeResults[,3:jmax]
	  # species names
	  nds <- names(nsr)	
	  nrs <- summary(nsr )
	  
	  # time info for nodes 
	  nodeTimeRangeSecs <- range(object$nodeResults$timeInSeconds)
	  nodeDeltaT <- mean(diff( object$nodeResults$timeInSeconds) )
	  
  } else {
	  # no node results 
	  nodeTimeRangeSecs <- NULL
	  nodeDeltaT <- NULL
	  nrs <- NULL 
	  nds <- NULL 
	  
  }

  ###############
  # link results
  ###############

  uli <- unique(object$linkResults$ID)
  if( length(uli) > 0 ){
	  lrn <- names(object$linkResults) 
	  jmax <- length(lrn) - 1 
	  lsr <- object$linkResults[,3:jmax]
	  lks <- names(lsr)
	  lrs <- summary(lsr )
	  linkTimeRangeSecs <- range(object$linkResults$timeInSeconds)
	  linkDeltaT <- mean(diff( object$linkResults$timeInSeconds) )
  } else {
	  # no links 
	  linkTimeRangeSecs <- NULL
	  linkDeltaT <- NULL 
	  lrs <- NULL 
	  lks <- NULL 
  }
  
  
  # collect into an object 
  msxrptSmry <- list( Title = object$Title,
                      #nodes 
                      numNodes = length( uni ), 
                      uniqueNodeIDs = uni,
					  nodeSpecies = nds,
                      nodeTimeRangeInSeconds = nodeTimeRangeSecs,
                      nodeTimestep = nodeDeltaT,
                      nodeResSmry = nrs,
                      #links
                      numLinks = length( uli ),
                      uniqueLinkIDs = uli,
					  linkSpecies = lks,
                      linkTimeRangeInSeconds = linkTimeRangeSecs,
                      linkTimestep = linkDeltaT,
                      linkResSmry = lrs ) 
  
  class(msxrptSmry) <- "summary.epanetmsx.rpt" 

  return( msxrptSmry)                     

}



#' Print msx rpt summary
#' 
#' The function prints a summary of multi-species simulation results
#' contained in the report file 
#'  
#' @export
#' @param x a summary.epanetmsx.rpt object   
#' @param ... further arguments passed to print 
print.summary.epanetmsx.rpt <- function(x,...){

  cat( x$Title)
  cat("\n") 

  # node results 
  cat(" node results\n")
  print( x$nodeResSmry)

  # link results 
  cat(" Link results\n")
  print( x$linkResSmry)

}



#' Plot method for epanetmsx.rpt
#' 
#' Plots a sparkline table of Epanet-msx results 
#' 
#' @export 
#' @param x epanetmsx.rpt object
#' @param elementType character indicating whether results for "nodes" or links" should be plotted 
#' @param ... further arguments passed to plotSparklineTable 
#' @seealso plotSparklineTable 
plot.epanetmsx.rpt <- function(x, elementType = 'Nodes',...){
	
	# argument checking 
	if( length(elementType) > 1 ) stop("elementType must have length 1")
	
	sx <- summary(x)
	
	if( grepl("nodes", elementType, ignore.case = TRUE ) ){
	
	   xl <- c(   utils::head(x$nodeResults$Time, 1),
			      utils::tail(x$nodeResults$Time, 1) ) 
		
		  plotSparklineTable( x$nodeResults, row.var = 'ID', col.vars = sx$nodeSpecies, 
				xvar = 'timeInSeconds',
				xrange.labels = xl , ... )
		
		
	} else if( grepl("links", elementType, ignore.case = TRUE) ){
		   xl <- c(   utils::head(x$linkResults$Time, 1),
			      utils::tail(x$linkResults$Time, 1) ) 
	
		plotSparklineTable( x$linkResults, row.var = 'ID', col.vars = sx$linkSpeices, 
				xvar = 'timeInSeconds',
				xrange.labels = xl, ... )
		
	} else { 
		
		stop("illegal value of argument 'elementType' ")
		
	}
}





 





