NULL

#' 
#'\code{[} S3 method for \code{\link{smet-class}} object

#' @param x a \code{\link{smet-class}} object 
#' @param rows rows of the \code{x} data which will be returned. 
#' @param fields column number or fields of \code{x} which will be returned.
#' @param file file of the new returened \code{\link{smet-class}}. Default is \code{NA}.
#' @param date.field field name of the Date and Time variable. Default is \code{"timestamp"}. If it is not \code{NA}, it is always returned.
#' @param ... further argument for \code{\link{[}} method
#' 
#' @export
#' @rdname extract
#' @method [ smet

#' @aliases [ Extract
# @usage \method{[}{blockmatrix} (M, i = 1:nrow(M), j = 1:ncol(M),numeric_value=TRUE,blockmatrix=FALSE,...) 
#' @author Emanuele Cordano 
#' 
#' @examples 
#' 
#' x <- smet(system.file('examples/PIEM001114.smet',package="RSMET"))
#' x[,1:3]
#' 
#' 





'[.smet' <- function (x,rows="all",fields="all",file=NA,date.field="timestamp",...)  {
	
	header_names <- c("fields","units_offset","units_multiplier")
	
	cond <- all(header_names %in% names(x@header) )
	
	if (is.null(date.field)) date.field <- NA
	
	if (cond==FALSE) {
		
		id <- x@header$station_id[1]
		msg <- sprintf("Malformed SMET object: %s",id)
		stop(msg) 
		
	}
	
	if (is.numeric(fields)) fields <- x@header$fields[fields]
	
	if (any(fields=="all")) fields <- x@header$fields

	cond <- all(fields %in% x@header$fields)
	
	if (cond==FALSE) {
		
		warning("Some fields are ignored bacause not present!")
		
		
		
	}
	
	fields <- unique(c(date.field,fields))
	fields <- fields[fields %in% x@header$fields]
	
	out <- x
	
	
	out@header$fields <- fields
	out@header$units_offset <- out@header$units_offset[fields]
	out@header$units_multiplier <- out@header$units_multiplier[fields]
	
	###
	
	if (rows[1]=="all") rows <- 1:nrow(out@data)
	out@data <- out@data[rows,fields]
	out@file <- as.character(file)
	###
	return(out)
	
}