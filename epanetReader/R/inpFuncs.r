#************************************
#
#  (C) Copyright IBM Corp. 2015
#
#  Author: Bradley J Eck
#
#************************************

#    File:  inpFuncs.r
#
#    By:  bradley.eck@ie.ibm.com
#
#    Purpose:  functions for reading .inp files into R 
#
     




.lineRange <- function(tag, allLines){
  # get the range of line numbers that pertain
  # to a table 
  #
  begin <- grep(tag, allLines)  + 1 
  
  if( length(begin) > 1 ){
    stop( paste(tag, "appeared more than once"))
  }
  # file lines starting with [ aka taglines  
  tl <- grep("\\s*\\[",allLines)
  #choose the first one after begin 
  end <- tl[ which(tl>begin)[1] ] - 1 

  return( c(begin,end) )
}



## Remove Comments and excess whitespace
##
## Processes lines from an .inp file by (a) removing 
## characters after semi-colon (;) and (b) removing
## leading and trailing whitespace  
##
## @param someLines character vector that probably corresponds to 
##                  a section of an .inp file 
## @return character vector 
## @author Bradley J. Eck 
.processCommentsAndClean <- function(someLines){  
  
  #  semi-colon, then any character (.) 
  # any number of times (*) 
  # until the end of the string ($).
  nocmt <- gsub( ";.*$", "", someLines)
  
  # remove leading and trailing whitespace
  nolead <- gsub("^\\s+","", nocmt )
  notail <- gsub("\\s+$","", nolead)
  
  # change remaining white space to a single space
  # so that char vects are more readable 
  clean <- gsub("\\s+", " ", notail )
  
  # just get lines that have something 
  goodLines = grep(".", clean)
  
  res <-  clean[goodLines ]  
  
  return( res )
}


.inpSection2char <- function(tag, allLines)
{
  #private helper function to extract a section
  # of an inp file, remove comments and empty lines
  # and store as character vector
  
  
  # range of lines for the tag of interest
  rg <- .lineRange(tag, allLines)
  
  if( is.na(rg[1]) ){ 
    # the section is not in the .inp file
    return( NULL )
  } else {
    # the section is there proceed as usual

    # get the lines that pertain to the table 
    someLines <- allLines[rg[1]:rg[2]]
    
    # prep by removing comments and empty lines
    preppedLines <- .processCommentsAndClean(someLines)
   
    # catch if these are empty 
    if( length(preppedLines) == 0 ){
      preppedLines <- NULL
    }

    return( preppedLines )
  }
}

.inpSection2df <- function(tag, allLines){
  # private helper function to get some section from
  # an inp file and read it to a data frame 
  
  sect <- .inpSection2char(tag, allLines)

  if( is.null(sect)){
    return(NULL)
  } else { 
  
    # convert the data into a data frame 
    df <- utils::read.table( text= sect, as.is = TRUE, 
                      fill = TRUE, header = FALSE)  
    
    return( df )
  }
}



## Junctions Table 
##
## Convert the [JUNCTIONS] section of an .inp file to a data frame 
##
## @param allLines output of readLines for the .inp file 
## @return A data frame with column names corresponding to those
##        exported from the Epanet GUI
JUNCTIONS <- function( allLines){
	
	tag <- "\\[JUNCTIONS\\]"
	df <- .inpSection2df(tag, allLines)  
	
	if( is.null( df ) ){
		return( NULL )  
	} else {   
		# now process the column names according to what was input
		names(df)[1] <- "ID"
		df$ID <- as.character(df$ID)
		
		names(df)[2] <- "Elevation"
		
		if( (dim(df)[2]) > 2 ){
			# there is a demand column already
			names(df)[3] = "Demand"
		} else {
			# add it anyway and fill with NA 
			df$Demand <- NA
		}
		
		if( dim(df)[2] > 3){
			# there is a pattern column
			names(df)[4] <- "Pattern"
			df$Pattern <- as.factor(df$Pattern)
		} else {
			# add it anyway and fill w NA
			df$Pattern <- NA
		}
		
		return(df)
	} 
}

## Reservoirs Table 
##
## Convert the [RESERVOIRS] section of an .inp file to a data frame 
##
## @param allLines output of readLines for the .inp file   
## @return A data frame with column names corresponding to those
##        exported from the Epanet GUI
RESERVOIRS <- function(allLines){
	
	tag <- "\\[RESERVOIRS\\]"
	df <- .inpSection2df(tag, allLines)  
	
	if( is.null( df ) ){
		return( NULL )  
	} else {   
		# now process the column names according to what was input
		names(df)[1] <- "ID"
		df$ID <- as.character(df$ID)
		
		names(df)[2] <- "Head"
		
		if( dim(df)[2] > 2){
			# there is a pattern column
			names(df)[3] <- "Pattern"
			df$Pattern <- as.factor(df$Pattern)
		} else {
			# add it anyway and fill w NA
			df$Pattern <- NA
		}
		
		return( df )
	}
	
}

## Tanks Table 
##
## Convert the [TANKS] section of an .inp file to a data frame 
##
## @param allLines name of the .inp file to read from 
## @return A data frame with column names corresponding to those
##        exported from the Epanet GUI
TANKS <- function( allLines ){
	tag <- "\\[TANKS\\]"
	df <- .inpSection2df(tag, allLines)  
	
	if( is.null(df)) {
		return( NULL )  
	} else {   
		# rename the columns
		names(df)[1:6] <- c("ID","Elevation","InitLevel",
				"MinLevel", "MaxLevel","Diameter")  
		
		#convert id field to character
		df$ID <- as.character(df$ID)
		
		# deal w optional fields
		if( dim(df)[2]>6){
			names(df)[7] <- "MinVol"        
		} else {
			df$MinVol = NA
		}
		
		if( dim(df)[2]>7 ){
			names(df)[8] <- "VolCurve"
			df$VolCurve <- as.factor(df$VolCurve)
		} else {
			df$VolCurve <- NA
		}
		
		return(df)
	} 
}


## Pipes Table 
##
## Convert the [PIPES] section of an .inp file to a data frame 
##
## @param allLines name of the .inp file to read from 
## @return A data frame with column names corresponding to those
##        exported from the Epanet GUI
PIPES <- function( allLines ){
	tag <- "\\[PIPES\\]"
	df <- .inpSection2df(tag, allLines)  
	
	
	if( is.null(df ) ) { 
		return( NULL  )  
	} else {   
		
		# rename the columns
		names(df)[1:6] <- c("ID","Node1","Node2","Length",
				"Diameter","Roughness")
		
		#convert id field to character
		df$ID <- as.character(df$ID)
		
		# Deal w optional fields
		if( dim(df)[2]>6){
			names(df)[7] <- "MinorLoss"        
		} else {
			df$MinorLoss = NA
		}
		
		if( dim(df)[2]>7 ){
			names(df)[8] <- "Status"
			df$Status <- as.factor(df$Status)
		} else {
			df$Status<- NA
		}
		
		return(df)
	}
}

## Pumps Table 
##
## Convert the [PUMPS] section of an .inp file to a data frame 
##
## @param allLines name of the .inp file to read from 
## @return A data frame with column names corresponding to those
##        exported from the Epanet GUI
PUMPS <- function( allLines ){
	tag <- "\\[PUMPS\\]"
	df <- .inpSection2df(tag, allLines)  
	
	
	if( is.null( df )) { 
		return( NULL )  
	} else {   
		# rename the columns
		names(df)[1:3] <- c("ID","Node1","Node2")
		
		#convert id field to character
		df$ID <- as.character(df$ID)
		
		#collapse any fields beyond Parameters into Parameters
		nc <- dim(df)[2]
		
		#not quite sure how this works, but it does!
		df$Parameters <- do.call(paste, df[,4:nc])
		
		# just keep the four cols we like
		pmp <- df[,c("ID", "Node1", "Node2", "Parameters")]
		
		# store the parameters as a factor 
		pmp$Parameters <- as.factor(pmp$Parameters)
		
		return(pmp)
	}
}

## Valves Table 
##
## Convert the [VALVES] section of an .inp file to a data frame 
##
## @param allLines name of the .inp file to read from 
## @return A data frame with column names corresponding to those
##        exported from the Epanet GUI or NULL if the table is missing
VALVES <-function( allLines){
  
  tag <- "\\[VALVES\\]"
  df <- .inpSection2df(tag, allLines)  
  
#  if( is.na( unlist(df)[1] )){
  if( is.null(df) ) { 
    return( NULL )  
  } else {   
    #proceed as usual
    
    # rename the columns
    names(df)[1:3] <- c("ID","Node1","Node2")
    
    #convert id field to character
    df$Node1 <- as.character(df$Node1) 
    df$Node2 <- as.character(df$Node2) 

    # name further cols 
    names(df)[4:7] <- c("Diameter", "Type", "Setting", "MinorLoss")
    df$ID <- as.character(df$ID)
	df$Type <- as.factor(df$Type)
  
    return(df)
  }
}


## Demands Section
##
## read [DEMANDS] section of an .nip file to a data.frame 
##
## @param allLines 
## @return data.frame of demands 
DEMANDS <-function( allLines){
  tag <- "\\[DEMANDS\\]"
  df <- .inpSection2df(tag, allLines)  
  
  if( is.null(df) ) { 
    return( NULL )  
  } else {   
    #proceed as usual
    
    # rename the columns
    names(df)[1:3] <- c("Node","Demand","Pattern")
    
    #convert id and pattern field to character
    df$Node <- as.character(df$Node) 
	df$Pattern <- as.factor(df$Pattern)

    return(df)
  }
}


## Energy Parameters
##
## Convert the [ENERGY] section of an .inp file to a  character vector
##
## @param allLines results of readLines on .inp file  
## @return A character vector with an entry for each line of the section or 
##         NA if the section is missing. 
##         Comments and excess whitespace are removed. 
ENERGY <- function( allLines ){
  tag <- "\\[ENERGY\\]"
  sect <- .inpSection2char(tag,allLines)
  return(sect)  
}

## Time Parameters
##
## Convert the [TIMES] section of an .inp file to a  character vector
##
## @param allLines name of the .inp file to read from 
## @return A character vector with an entry for each line of the section or 
##         NA if the section is missing. 
##         Comments and excess whitespace are removed. 
TIMES <- function( allLines ){
  tag <- "\\[TIMES\\]"
  sect <- .inpSection2char(tag,allLines)
  return( sect )
}

## Options 
##
## Determine the options specified by the [OPTIONS] section of an .inp file
## @details The [OPTIONS] section is read from the allLines and used to update
##          a list of Epanet's default options. In this way if an option such as
##          units is not specified by the .inp file, the units that would be used by
##          default are provided.
##
## @param allLines name of the .inp file to read from 
## @return A list with named entries for the options of Epanet.
## @references Rossman, L. A. (2000). Epanet 2 users manual. US EPA, Cincinnati, Ohio.
## http://nepis.epa.gov/Adobe/PDF/P1007WWU.pdf
OPTIONS <- function( allLines ){
  
tag <- "\\[OPTIONS\\]"
  sect <- .inpSection2char(tag,allLines)
  
  # update with changes as specified in the file 
  updatedOptions <- .listUpdater( epanetDefaultOptions(), sect)
  
  return( updatedOptions ) 
}

#' Epanet Default Options
#' 
#'  A list of Epanet's default options
#'
#' @export 
#' @details Provides a named list in the form of OPTION = default_value where the 
#' values are taken from pages 152-154 of the manual.
#' @references Rossman, L. A. (2000). Epanet 2 users manual. US EPA, Cincinnati, Ohio.
#' http://nepis.epa.gov/Adobe/PDF/P1007WWU.pdf
#' @examples 
#' epanetDefaultOptions() 
epanetDefaultOptions <- function(){
  defaultOptions <- list(
    UNITS = "GPM",
    HEADLOSS ="H-W",
    HYDRAULICS = NA,
    QUALITY = "NONE",
    VISCOSITY = 1.0,
    DIFFUSIVITY = 1.0,
    SPECIFIC_GRAVITY = 1.0, 
    TRIALS = 40,
    ACCURACY = 0.001,
    UNBALANCED = "STOP",
    PATTERN = NA,
    DEMAND_MULTIPLIER = 1.0,
    EMITTER_EXPONENT =  0.5,
    TOLERANCE = 0.01,
    MAP = NA )
  
  return(defaultOptions)
}


## Update a list with a character vector
##
## Change values in a named list 
## based on info in a character vector.
## Assumes  multi-word names in charVec are joined
## by underscore in the named list.
.listUpdater <-function( oldList, charVec ){
  
  imax <- length( charVec ) 
  for( i in 1:imax ){
    
    tokens <-  unlist(strsplit(charVec[i], "\\s+"))
    nt <- length(tokens)  
    # last token is the value 
    val <- tokens[nt]
    # remaining tokens form the key / name 
    key <- paste(tokens[1:(nt-1)], collapse="_")
    
    # find the list entry pertaining to this key 
    list_entry <- grep( key, names(oldList), ignore.case=TRUE)
    lle <- length( list_entry) 
    
    # update the list if there is one possible entry
    if( lle == 1 ){
    
      #update the list
      oldList[[list_entry]] <- val
    } else if ( lle > 1 ){
      # throw an error if multiple options are possible 
      stop(paste("associated multiple options with", charVec[i]))
    }
    
  }
  
  return(oldList)  
}



## Coordinates Table 
##
## Convert the [COORDINATES] section of an .inp file to a data frame 
##
## @param allLines name of the .inp file to read from 
## @return A data frame with column names corresponding to those
##        exported from the Epanet GUI or NULL if the table is missing
COORDINATES <- function( allLines ){
  
  tag <- "\\[COORDINATES\\]"
  df <- .inpSection2df(tag, allLines)  

  if( is.null(unlist(df)[1])){
    return( NULL )  
  } else {   
    #proceed as usual
    
    # rename the columns
    names(df)[1:3] <- c("Node","X.coord","Y.coord")
    
    #convert id field to character
    df$Node <- as.character(df$Node)
    
    return( df )
  }
}

## Pattern Information   
##
## Convert the [PATTERNS] section of an .inp file to a named list
## @details Pattern information is extracted from the .inp file and 
##          returned as a list with entries named as the ID of the pattern.
##
##          If integers are used as pattern IDs names of list elements are backquoted
##          according to the default behavior in R.  So if the .inp file has a pattern "1"
##          this pattern will appear as element `1` in the list that is returned.
##
## @param allLines name of the .inp file to read from 
## @return A named list or NA if the section is missing 
PATTERNS <- function( allLines){
  
  tag <- "\\[PATTERNS\\]"
  sect <- .inpSection2char(tag,allLines)
 
  if( is.null(sect) ) {
    return (NULL)
  } else { 

    # break lines of section into tokens
    L <-  strsplit(sect, "\\s+")

    # the first token of each line is the pattern ID 
    firstToken <- unlist(lapply( L , function(x) x[1]))
    patternIDs <- unlist(unique(firstToken)  )

    # for each pattern ID, get the requisite lines of L

    .getPattern <- function(ID){

      # rows that pertain to this pattern
      rows <-  which( firstToken == ID )

      # second thru last entry on these lines 
      LL <-lapply( L[rows], function(x) x[2:length(x)] )  

      # combine these together 
      pat <- as.numeric(unlist( LL ) )

    }

    patterns <- lapply( patternIDs, .getPattern) 

    names(patterns) <- patternIDs    

    #warning if patternIDs contain integers
    anyints <- as.logical( max( grepl("^[0-9]", patternIDs ) ) )
    if( anyints == TRUE ){
       warning("patterns have integer IDs, see ?epanet.inp")
    }
   
    return(patterns) 
  }
  
}

## Curve Information   
##
## Convert the [CURVES] section of an .inp file to a named list
##
## @details Curve information is extracted from the .inp file and 
##          returned as a list with entries named as the ID of the curves.
##
##          If integers are used as curve IDs names of list elements are backquoted
##          according to the default behavior in R.  So if the .inp file has a curve "1"
##          this pattern will appear as element `1` in the list that is returned.
##
## @param allLines name of the .inp file to read from 
## @return A named list or NULL if the section is missing 
CURVES <- function(allLines){

  tag <- "\\[CURVES\\]"
  sect <- .inpSection2char(tag,allLines)

  if( is.null(sect) ) {
    return (NULL)
  } else { 
    # break lines of section into tokes
    L <-  strsplit(sect, "\\s+")

    # the first token of each line is the pattern ID 
    firstToken <- unlist(lapply( L , function(x) x[1]))
    IDs <- unlist(unique(firstToken)  )

    # for each ID, get the requisite lines of L
    .getCurve <- function(ID){

      # rows that pertain to this pattern
      rows <-  which( firstToken == ID )

      # X values   
      X <-unlist(lapply( L[rows], function(x) x[2] )  )
      # Y values   
      Y <-unlist(lapply( L[rows], function(x) x[3] )  )


      # combine these together 
      curve <- list(X=as.numeric(X),Y=as.numeric(Y))

      return(curve)

    }

    curves <- lapply( IDs, .getCurve) 
    names(curves) <- IDs

    #warning if curve IDs contain integers
    anyints <- as.logical( max( grepl("^[0-9]", IDs ) ) )
    if( anyints == TRUE ){
       warning("curves have integer IDs, see ?epanet.inp")
    }
    return(curves)
  }
}

## Title
## 
## Read the [TITLE] section of an .inp File
## 
## @param allLines  path of the .inp file to read in 
TITLE <- function( allLines ){
  
  tag <- "\\[TITLE\\]"
  sect <- .inpSection2char(tag,allLines)
  return( sect )
}


STATUS <- function( allLines ){
	
	tag <- "\\[STATUS]"
    df <- .inpSection2df(tag, allLines )
	
	if( is.null( df ) ){
		return( NULL )  
	} else {   
		# now process the column names according to what was input
		names(df)[1] <- "ID"
		df$ID <- as.character(df$ID)
		
		names(df)[2] <- "Status"
        df$Status <- as.factor(df$Status)

		return(df)
	}
}


