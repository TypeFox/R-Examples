
NULL

#' Function to read a  a SMET file
#' 
#' 
#' 
#' 
#' @param file  SMET file name
#' @param numeric logical value
#' @param timezone.offset.sign It can be \code{"negative"} or \code{"positive"}. Many systems support timezones of the form GMT+n and GMT-n, which are at a fixed offset from UTC (hence no DST). Contrary to some usage (but consistent with names such as PST8PDT), negative offsets are times ahead of (east of) UTC, positive offsets are times behind (west of) UTC.
#' @param date.field field neme used for date and time. Default is \code{"timestamp"}.
#' @param date.format format used for date and time. Default is \code{"\%Y-\%m-\%dT\%H:\%M:\%S"}.
#' @param non_numeric_fields fields, except \code{date.fiels}, that contain non-numeric values
#' @param comment character symbol used to comment lines or part of lines. Default is \code{"#"}. See \code{\link{str_locate}}
#' @param station.field field name used for station ID. Default is \code{"station_id"}, as used for \code{SMET} format.
#' @param ... further arguments
#' 
#' @export
#' @import stringr
#' 
#' @return a \code{\link{smet-class}} object
#' 
#' @details To better understand the use of timezones in R (\code{tz} attribute) see the following link: 
#' \url{http://stackoverflow.com/questions/11927433/timezones-in-r-how-to-avoid-ambiguous-terms-such-as-est}
#' 
#' @seealso \code{\link{smet-class}}
#' 
#' @examples 
#' 
#' file <- system.file("examples/test.smet",package="RSMET")
#' sm <- smet(file)
#' 

smet <- function(file=NULL,numeric=TRUE,non_numeric_fields=NULL,
		timezone.offset.sign=c("negative","positive","-","+"),
		date.field="timestamp",
		date.format="%Y-%m-%dT%H:%M:%S",comment="#",station.field="station_id",
		...) {
	
	
	
	
	file_default <- system.file("examples/test.smet",package="RSMET")
	
	
	
	
	if (is.null(file)) file <- NA
	if (is.na(file)){
		
		
		file <- file_default
		warning("file is missing, and automatically set by default!")
		
	}
	
	
	
	string <- readLines(file,encoding="US-ASCII")
	
	#####
	
	sl <- stringr::str_length(string)
	sla <- stringr::str_locate(string,comment)[,"start"]-1
	sla[is.na(sla)] <- sl[is.na(sla)]
	string <- stringr::str_sub(string,0,sla)
	string <- string[string!=""]
	
	#####
	
	
	
	iheader <- which(string=="[HEADER]")
	idata <- which(string=="[DATA]")
	
	signature <- paste(string[1:(iheader-1)],collapse=";")
	
	outdata <- stringr::str_split(string[-(1:idata)]," ")
	outdata <- lapply(X=outdata,FUN=function(x){ x[x!=""]}) ## ec 20151118
	
	outdata <- t(as.data.frame(outdata,stringsAsFactors=FALSE))
	data <- as.data.frame(outdata,stringsAsFactors=FALSE)
	rownames(data) <- NULL
	
	#### Create a list for the header
    outheader <- stringr::str_split(string[(iheader+1):(idata-1)],"=",n=2)                                                             
	header <- lapply(X=outheader,FUN=function(x){x[2]})
	names(header) <- sapply(X=outheader,FUN=function(x){x[1]})
	names(header) <- stringr::str_replace_all(names(header)," ","")
	
	
	#### Data and Header Processing
    numerickeys <- c("nodata","altitude","latitude","longitude","easting","northing")	

	numerickeys <- numerickeys[numerickeys %in% names(header)]
	for (it in numerickeys) {
		
		header[[it]] <- as.numeric(header[[it]])
	}
	
	vars <- header$fields
	vars <- (stringr::str_split(vars," ")[[1]])
	vars <- vars[vars!=""]
	header$fields <- vars
	
	
	if ("units_offset" %in% names(header)) {
		
		offset <- header[["units_offset"]]
		
		offset <- (stringr::str_split(offset," ")[[1]])
		
		offset <- as.numeric(offset[offset!=""])
	
		header[["units_offset"]] <- offset
		names(header[["units_offset"]]) <- header$fields
		
	} else {
		
		header[["units_offset"]] <- array(0,ncol(data))
		names(header[["units_offset"]]) <- header$fields
		
		
	}
	
	if ("units_multiplier" %in% names(header)) {
		
		mult <- header[["units_multiplier"]]
		
		mult <- (stringr::str_split(mult," ")[[1]])
		
		mult <- as.numeric(mult[mult!=""])
		header[["units_multiplier"]] <- mult
		names(header[["units_multiplier"]]) <- header$fields
		
	} else {
		
		header[["units_multiplier"]] <- array(1,ncol(data))
		names(header[["units_multiplier"]]) <- header$fields
		
		
	}
	
	names(data) <- header$fields
	names(header$fields) <- header$fields
	#### END Data and Header Processing
	
	if (numeric==TRUE) {
		
		stringsAsFactors=FALSE
		non_numeric_fields <- c(non_numeric_fields,date.field)
		
		numeric_fields <- names(data)[!(names(data) %in% non_numeric_fields)]
		
		for (it in numeric_fields) {
			
			data[,it] <- as.numeric(data[,it])
		}
		
		timezone.offset.sign <- timezone.offset.sign[1]
	
		if (timezone.offset.sign %in% c("positive","+")) timezone.offset.sign <- 1
		if (timezone.offset.sign %in% c("negative","-")) timezone.offset.sign <- -1
		
		if (timezone.offset.sign!=0) timezone.offset.sign <- timezone.offset.sign/abs(timezone.offset.sign)
		
		tz <- header$tz
		
		if (is.null(tz)) {
			
			tz <- base::list(...)$tz
		} 
		if (is.null(tz)) {
			tz <- NA
			
		}
		
		tzn <- as.numeric(tz)*timezone.offset.sign		
		
		
		if (!is.na(tzn)) {
			
			
			if (tzn>=0) {
				
				tz <- sprintf("Etc/GMT+%d",tzn)
				
				
			} else {
				
				
				tz <- sprintf("Etc/GMT-%d",-tzn)
				
			}
			
			
		}
		
		if (date.field %in% names(data)) {
			
			data <- data
			times <- data[,date.field]
			data <- data[,names(data)!=date.field]
			
			
			data$times_temp <- as.POSIXlt(times,format=date.format,tz=tz)
			
			
			names(data)[names(data)=="times_temp"] <- date.field
			data <- data[,header$fields]
			
			####
			
	#		> (df$timestamp)
	#		[1] "2015-10-21 07:20:00 GMT" "2015-10-22 07:30:00 GMT"
	#		[3] "2015-10-22 11:00:00 GMT" "2015-10-23 11:00:00 GMT"
	#		[5] "2015-10-24 07:30:00 GMT" "2015-10-24 11:00:00 GMT"
	#		[7] "2015-10-25 07:30:00 GMT" "2015-10-27 08:00:00 GMT"
	#		[9] "2015-10-27 12:00:00 GMT" "2015-10-28 07:40:00 GMT"
	#		[11] "2015-10-29 12:00:00 GMT" "2015-10-29 08:30:00 GMT"
	#		[13] "2015-10-31 08:00:00 GMT" "2015-11-03 08:30:00 GMT"
	#		[15] "2015-11-04 08:30:00 GMT" "2015-11-04 12:00:00 GMT"
	#		> sort(df$timestamp)
	#		[1] "2015-10-21 07:20:00 GMT" "2015-10-22 07:30:00 GMT"
	#		[3] "2015-10-22 11:00:00 GMT" "2015-10-23 11:00:00 GMT"
	#		[5] "2015-10-24 07:30:00 GMT" "2015-10-24 11:00:00 GMT"
	#		[7] "2015-10-25 07:30:00 GMT" "2015-10-27 08:00:00 GMT"
	#		[9] "2015-10-27 12:00:00 GMT" "2015-10-28 07:40:00 GMT"
	#		[11] "2015-10-29 08:30:00 GMT" "2015-10-29 12:00:00 GMT"
	#		[13] "2015-10-31 08:00:00 GMT" "2015-11-03 08:30:00 GMT"
	#		[15] "2015-11-04 08:30:00 GMT" "2015-11-04 12:00:00 GMT"
	#		> help(sort,help_type="html")
			row.names(data) <- sprintf("%09d",1:nrow(data))
			times <- data[,date.field]
			names(times) <- row.names(data)
			
			data <- data[names(sort(times)),]
			
						
			####
			
			
			
			
		}
	}
	
	if (is.null(file)) file <- NA
	if (identical(file,file_default)) file <- NA
	
	if (station.field %in% names(header)) {
		
		station_id <- header[[station.field]][1]
		
		station_id <- str_replace(station_id," ","")
		
		header[[station.field]] <- station_id
	}
	
	out <- new("smet",signature=signature,header=header,data=data,file=as.character(file))
	
	
	
	
	
	return(out)
	
	
	
	
	
}