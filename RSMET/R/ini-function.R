
NULL

#' Function to read a  a Meteo IO *.ini file
#' 
#' @param file  Meteo IO INI file name
#' @param comment character symbol used to comment lines or part of lines. Default is \code{"#"}. See \code{\link{str_locate}}
#' @param ... further arguments
#' 
#' @export
#' @import stringr
#' 
#' @return a \code{\link{meteoioini-class}} object
#' 

#' @seealso \code{\link{meteoioini-class}},\code{\link{print.meteoioini}}
#' 
#' @examples 
#' 
#' file <- system.file("examples/io.ini",package="RSMET")
#' ini <- meteoioini(file)
#' 

meteoioini <- function(file=NULL,comment="#",
		...) {
	
	
	
	object_name <- "meteoioini"
	file_default <- system.file("examples/io.ini",package="RSMET")
	
	out <- new(object_name)
	
	
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
	
	
	
	slots <- getSlots(object_name)
	
	inienvs <- names(slots[slots=="list"])
	inienvsbr <- paste("[",inienvs,"]",sep="")
	
	################
	names(inienvs) <- inienvs
	inienvs <- inienvs
	ub <- unlist(sapply(X=inienvs,FUN=function(x,string) { which(stringr::str_detect(string,x))},string=string))
###	ubn <- inienv[ub]
    ii <- findInterval(1:length(string),ub)
	groups <- names(ub)[ii]
	
	
##	ma <- cbind(string,findInterval(1:length(ini),ub)])
	
	for (it in unique(groups)) {
		
		##print(it)
		value <- string[groups==it]
		value <- value[stringr::str_detect(value,"=")]
		value <- as.list(value)
		value <- stringr::str_split(value,"=")
		names(value) <- lapply(X=value,FUN=function(x){x[1]})
		value <-  lapply(X=value,FUN=function(x){x[-1]})
	
		names_v <- str_replace_all(names(value)," ","")
		names_v <- str_replace_all(names_v,"[\t]","") ## 		https://www.google.it/webhp?sourceid=chrome-instant&ion=1&espv=2&ie=UTF-8#q=%5Ct
		names_v <- str_replace_all(names_v,"[\b]","") ## 		https://www.google.it/webhp?sourceid=chrome-instant&ion=1&espv=2&ie=UTF-8#q=%5Ct
	
		
		names(value) <- names_v
		
		value <- lapply(X=value,FUN=stringr::str_replace_all,pattern="[\t]",replacement=" ")
		value <- lapply(X=value,FUN=stringr::str_replace_all,pattern="[\b]",replacement=" ")
		value <- lapply(X=value,FUN= stringr::str_split,pattern=" ")
		value <- lapply(X=value,FUN=function(x){x[[1]]})
		value <- lapply(X=value,FUN=function(x){x[x!=""]})
		
		slot(out,it) <- value
		
	}

	
	if (identical(file,file_default)) {
		
		out@file <- as.character(NA)
	} else {
		
		out@file <- file
	}
	
	return(out)
	
	
	
}