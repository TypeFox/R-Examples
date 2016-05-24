NULL

#' 
#'\code{'[<-'} S3  Replacement method for \code{smet-class} object
#' 
#' @param x a \code{\link{smet-class}} object 
#' @param rows rows of the \code{x} data which will be replaced. 
#' @param fields column number or fields of \code{x} which will be replaced.
#' @param nodata numeric value used to replce \code{NA}. It is used if no-data value is not specified in the header of \code{x}.
#' @param date.field field name of the Date and Time variable. Default is \code{"timestamp"}. If it is not \code{NA}, it is always returned.
#' @param multiplier scalar value of the unit multiplier for each new \code{field} name which is not previously present.
#' @param offset scalar value of the unit multiplier for each new \code{field} name which is not previously present.
#' @param ... further argument for \code{\link{[}} method
#' 
#' 
#' 
#' 
#' @param value a object to be replaced
#' @export
#' 
#' @return The "replaced" \code{\link{smet-class}} object.
#' @note In case \code{fields} is a character vector, the elements whose names is in \code{value} is replaced.
#' @rdname extract_replacemethod
#' @method [<- smet

#' @aliases [<-,extract_replacemethod

#' @author Emanuele Cordano 
#' 
#' @examples 
#' rm(list=ls())
#' 
#' 
#'  ### SMET Rifugio Vaccarone
#' @examples 
#' 
#' x <- smet(system.file('examples/PIEM001114.smet',package="RSMET"))
#' x[,1:3]
#' 
#' x[,"VAR07"] <- NA
#' x[,"VAR08"] <- 123
#'  x
#' 
#' y <- x 
#' 
#' y[,1:2] <- x[,1:2]
#' 
#' y[,c("TSS1","TSS2"),offset=273.15] <- 0
#' 
 







'[<-.smet' <- function (x,rows="all",fields="all",date.field="timestamp",offset=0,multiplier=1,nodata=-9999,...,value)  {
	
	
	if (class(value)=="smet") {
		
		
		multiplier <- units_multiplier(value)
		offset     <- units_multiplier(value)
		value      <- as.data.frame(value)
		
		
		
	}
	
	
	
	header_names <- c("fields","units_offset","units_multiplier")
	
	cond <- (header_names %in% names(x@header) )
	
	if (is.null(date.field)) date.field <- NA
	
	if (any(cond==FALSE)) {
		
		
		id <- x@header$station_id[1]
		mis <- paste(cond[cond==FALSE],collapse=",")
		msg <- sprintf("Malformed SMET object: %s (%s missing)",id,mis)
		stop(msg) 
		
	}
	
	if (is.numeric(fields)) fields <- x@header$fields[fields]
	if (any(fields=="all")) fields <- x@header$fields
	
	oldf <-(fields %in% x@header$fields)
	##
	multiplier <- array(multiplier,length(fields))
	offset <- array(offset,length(fields))
	names(multiplier) <- fields
	names(offset) <- fields
	names(fields) <- fields
	
	
	multiplier[fields[oldf]] <- x@header$units_multiplier[fields[oldf]]
	offset[fields[oldf]] <-  x@header$units_offset[fields[oldf]]
	names(x@header$fields) <- x@header$fields
	
	
	
#	if (cond==TRUE) {
#		
#		warning("Some fields are ignored bacause not present!")
#		
#		added <- fiel
#		
#		
#		
#	}
	
#	fields <- unique(c(date.field,fields))
	######fields <- fields[fields %in% x@header$fields]
	
	out <- x
	
	
	out@header$fields[fields] <- fields
	
	
	out@header$units_offset[fields] <- offset[fields]
	out@header$units_multiplier[fields] <- multiplier[fields]
	
	
	
	###
	
	if (any(rows=="all")) {
		
		
		
		
		
		out@data[,fields] <- value
	} else { 
		
		
	
		out@data[rows,fields] <- value
	
	}	
		
	#####out@file <- as.character(file)
	###
	
	if (!is.null(out@header$nodata)) {
		
		nodata <- out@header$nodata[1]
		
		
	}
	
	out@data[is.na(out@data)] <- nodata
	
	return(out)
	
	
}	
