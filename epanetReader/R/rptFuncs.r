#************************************
#
#  (C) Copyright IBM Corp. 2015
#
#  Author: Bradley J Eck
#
#************************************

#  File: rptFuncs.r
#
#  By: bradley.eck@ie.ibm.com
#
#  Purpose: functions to read an rpt file 

.getSectionRange <- function( j, resLines, lengthOfAllLines ) {
  # get the section of the rpt file that corresponds 
  # to the j-th entry in resLines 
   
  startLine <- resLines[j]
  
  numTables <- length( resLines ) 
  
  if( j < numTables ){
    endLine <- resLines[j+1] - 1 
  } else {
    endLine <- lengthOfAllLines - 1 
  }
  
  rl <- list( start = startLine, end = endLine )  
  return( rl ) 
}

.getTimeStamp <- function( tsLine){
  
 tokens <- unlist(  strsplit(tsLine, "\\s") ) 

 stamp <-
   tokens[ grep(":.{2}:", tokens) ]

 if( length(stamp) > 1 ){
   msg <- paste( "more than one timestamp found on line", tsLine)
   stop( msg )
 }
 
 if( length(stamp) == 0 ){
	 # use a time of zero if none is supplied 
	 stamp <- "0:00:00"
 } 
 return( stamp ) 
 
}

.timeStampToSeconds <- function( stamp ){

	tokens <- unlist(strsplit(stamp,":"))
	
	if( length(tokens) < 2 ){
		stop(paste("don't have interpretation for stamp", stamp ))
	}
	if( length(tokens)  > 3 ){
		stop(paste("don't have interpretation for stamp", stamp ))
	}
	
	
	hours <- tokens[1]
	minutes <- tokens[2]

	m <- as.integer(minutes)
	if( m > 59){
		stop(paste0("invalid minutes in ", stamp))
	}
	
	
	totalSeconds <- as.integer(hours) * 3600 + 
			        m * 60 
	
	if( length( tokens)  == 3 ){ 
		
		seconds <- tokens[3]
		s <- as.integer(seconds)
		if( s > 59){
			stop(paste0("invalid seconds in ", stamp))
		}
		
		totalSeconds <- totalSeconds +  as.integer(seconds) 
	}
	
	return( as.integer(totalSeconds) ) 
}

.section2df <- function( sect ){
 # make a section into it's own data frame

  imax <- length(sect)
  
   # Take column headings from labels. Sometimes first column is labeled
   # in two rows, sometimes in 1.  See tests. 
   headerRow1 <- unlist(strsplit( gsub("^\\s+", "", sect[3] ), "\\s+"))
   lh1 <- length( headerRow1 )
   headerRow2 <- unlist(strsplit( gsub("^\\s+", "", sect[4] ), "\\s+"))
   lh2 <- length( headerRow2) 


   if( lh1 ==  lh2 ){
     columnNames <- c( headerRow1, "note" ) 
   } else if( lh2 == ( lh1 + 1 ) ) {
     columnNames <- c( headerRow2[1], headerRow1, "note" )
   } else {
     warning("unexpected header format in rpt file, check results")
     columnNames <- c( headerRow2[1], headerRow1, "note" )
   }

   # name the first column "ID" rather than
   # "Node" or "Link"  to be consistent with inp objects
	columnNames[1] <- "ID"
   
  # set colClasses, everything is numeric execpt fist 
	# and last column
	lcn <- length( columnNames)
  cc <- rep("numeric", lcn )
  cc[1] <- "character"  # for ID column 
  cc[lcn] <- "character" # for note column 
  
   
  # make the section a data frame 
  df <- utils::read.table( text = sect[6:imax], 
		  col.names=columnNames,
		  colClasses = cc, 
		  strip.white = TRUE, fill = TRUE, header = FALSE )
  
  # now add the time info to the table 
  stamp <- .getTimeStamp( sect[1] )  
  
  # convert that stamp into sections 
  time_secs <-  .timeStampToSeconds(stamp)
  
  df$Timestamp <- stamp
  df$timeInSeconds <- time_secs
  
  return( df ) 
}

#' Bin Breaker    
#' 
#' Generate break points for use with cut() 
#' and range labels based on sample max and min  
#' 
#' @param x vector to find cuts for 
#' @param nbin number of bins 
#' @return list with elements Breaks and Labels
#' @details 
#' Helpful in making labels use the acutal max and min rather than
#' the +/- 1% cut() uses by default. 
binBreaker <- function( x, nbin){
	
    xmax <- max(x, na.rm = TRUE )
    xmin <- min(x, na.rm = TRUE )
	
	# use cut() to find break points for the specified number of bins   
	labs <- levels( cut( x, breaks = nbin ))
	# from help(cut) example on getting break points 
	brks <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
				  upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
    
    # get the breaks to use with cut()   	
	brkpts4cut <-  c( brks[1,1], brks[,2])
	
	# force a zero in the legend values  	  
	labs4legend <- levels( cut( x, breaks = c( xmin, brks[1:(nbin-1),2], xmax) ))  
	
	#combine results into list and return 
	return( list( Breaks = brkpts4cut, Labels = labs4legend))
		
} 



checkRptFile <- function( allLines ){
  # check rpt file format 
   
   # look for new page character in the file 
   hasPageBreaks <-  as.logical(   max( grepl("\f", allLines) ) )

   if( hasPageBreaks ) {
      msg <- paste( " Page breaks not allowed in rpt file.\n",
                    "Put the line 'Page 0' in the [REPORT] section of\n",
                    "the .inp file and generate the .rpt file again.\n")
     stop( msg ) 
   }


   # look for node results and link results 
   hasNodeResults <- as.logical( max( grepl("Node Results", allLines)))
   hasLinkResults <- as.logical( max( grepl("Link Results", allLines)))

   if( hasNodeResults == FALSE ){
  
     msg <- paste(" Node results not found in .rpt file. \n",
                  "Add the line 'Nodes All' to the [REPORT] section of the .inp file.")
     warning(msg)
   }

   if( hasLinkResults == FALSE ){
  
     msg <- paste(" Link results not found in .rpt file. \n",
                  "Add the line 'Links All' to the [REPORT] section of the .inp file.")
     warning(msg)
   }


   if( ( hasNodeResults == FALSE )& 
       ( hasLinkResults == FALSE )  ){
       
       # no results to read, give error 
       stop("No results to read")

   }

}
