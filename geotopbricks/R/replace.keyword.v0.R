NULL

#'
#'It replaces some keyword values of \code{geotop.inpts} file with the ones of anoter \code{*.inpts}  value
#' 
#' @param x filename of the \code{*.inpts} with the "new" keyword value 
#' @param y filename of the \code{*.inpts} with the "old" keyword value. Default is \code{"geotop.inpts"}.
#' @param file.output filename where to write the comprehensive new \code{geotop.inpts} file. If it is \code{NULL} (default), the fileneme is assigned by \code{y}.
#' @param write.file.output logical value. If it is \code{TRUE}, the output of the function is written in he file \code{file.output}.
#' @param wpath working path to the GEOtop simulation folder containing the \code{x} and \code{y} files.
#' @param ... further arguments
#' 
#' @author Emanuele Cordano
#' 
#' @details This function repleces some keword values of \code{y} with the ones indicated in \code{y}. It is useful to replace the meteo station metedata, for instance, when the meteorological station of a study cases are modified. 
#' The function returns the new \code{geotop.inpts} file as a vector of character strings. If \code{write.file.output==TRUE}, the output is written in an extarnal file, e.g. \code{"geotop.inpts"} newly (this option is suggested).
#' 
#' 
#' @export
#' 
#' @examples
#' 
#' library(geotopbricks)
#' wpath <- system.file('template/meteo_ex',package="geotopbricks")
#' x <- "meteo.inpts"
#' zl <- replace.keyword(x,wpath=wpath,write.file.output=FALSE)
#' 
#' 
#' 


replace.keyword <- function(x,y="geotop.inpts",file.output=NULL,write.file.output=TRUE,wpath=NULL,...) {
	
	if (!is.null(wpath)) {
		len <- str_length(wpath)
		
		if ((str_sub(x,len,len)=='/')) x <- str_sub(x,2)
		if ((str_sub(y,len,len)=='/')) y <- str_sub(y,2)
		
		x <- paste(wpath,x,sep="/")

		y <- paste(wpath,y,sep="/")


	}

	
	xl <- readLines(x,...)
	yl <- readLines(y,...)
	
	
	
	keywords <- declared.geotop.inpts.keywords(wpath=NULL,inpts.file=x,...)$Keyword	
	
	for (key in keywords) {
		
		 ikey <- str_detect(yl,key)
		 yl <- yl[!ikey]
		
	}
	
	
	zl <- c(yl,xl)
	
	if (is.null(file.output)) file.output <- y
	
	
	if (write.file.output) writeLines(zl,con=file.output,sep="\n")

		
	
	return(zl)
	
	
	
	
}