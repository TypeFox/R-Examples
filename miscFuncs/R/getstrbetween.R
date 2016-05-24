##' getstrbetween function
##'
##' A function used in web scraping. Used to simplify the searching of HTML strings for information.
##'
##' @param linedata a string 
##' @param start integer, where to start looking in linedata
##' @param startmark character string. a pattern identifying the start mark
##' @param endmark character string. a pattern identifying the end mark
##' @param include include the start and end marks?
##' @return the first string after start and between the start and end marks
##' @export

getstrbetween <- function(linedata,start,startmark,endmark,include=FALSE){ 	# returns first occurrence of string between "startmark" and "endmark" after "start" characters of the string 
	if (include){
		# includes start and end marks
		begin <- min(gregexpr(startmark,linedata)[[1]][gregexpr(startmark,linedata)[[1]]>start])
		end <- min(gregexpr(endmark,linedata)[[1]][gregexpr(endmark,linedata)[[1]]>begin])
	}
	else{
		# removes start and end marks
		begin <- min(gregexpr(startmark,linedata)[[1]][gregexpr(startmark,linedata)[[1]]>start]) + nchar(startmark)
		end <- min(gregexpr(endmark,linedata)[[1]][gregexpr(endmark,linedata)[[1]]>begin]) - 1
	}
	reqtext <- substr(linedata,begin,end)
	return(reqtext)
}
