# TODO: Add comment
# 
# Author: ecor
###############################################################################

NULL

#'
#' Imports a brick of raster ascii maps into a 'brick' object 
#' 
#' 
#' 
#' @param file fileneme of the 'brick' files containing the decimal formatter. It is \code{NULL} by default, otherwise it replaces  \code{file_suffix}, \code{formatter} and \code{file_extension}.
#' @param file_prefix character string suffix name of  the 'brick' files.
#' @param formatter string value.  Default is \code{"\%04d"} . 
#' @param file_extension strinf value. Default is \code{".asc"}
#' @param nlayers number of layers
#' @param use.read.raster.from.url logical value. Default is \code{FALSE}. (this is recommended in this function). If \code{TRUE} the RasterLayer are read with \code{\link{read.raster.from.url}}, istead of \code{\link{raster}} (otherwise). It is recomended in case the files whose paths are contained in \code{x} are remote and are 'http' addresses. In this cases the stand-alone method \code{raster(x)} does not always work and \code{use.read.raster.from.url} is necessary.  
#' @param start.from.zero logical value. Default is \code{FALSE}. If \code{TRUE} the formatter starts from \code{0000}, otherwise it starts from \code{0001}.  
#' @param crs coordinate system see \code{\link{RasterBrick-class}},\code{\link{brick}}, Default is \code{NULL}.
#' 
#' 
#' @export
#' 
#' @return the output is returned as a \code{\link{RasterBrick-class}} object
#'
#' 
#' @examples
#' library(geotopbricks)
#' library(raster)
#' file <- system.file("doc/examples/snowthickness",package="geotopbricks")
#' file <- paste(file,"SnowThickness0000L%04d.asc",sep="/")
#' # nlayers=15
#' nlayers <- 6 ## Only 6 layers are read to minimize the elapsed time of the example!!
#' b <- brick.decimal.formatter(file=file,nlayers=nlayers)
#' nlayers(b)
#' names(b)


# @examples
# map <- '~/LAVORI_IN_CORSO/2012/SNOWMAPS/openda_student/openda_student/exercise_geotop/stochModel/original_model/bin/bondone/rec/SnowThickness0004L%04d.asc'
#   
#




brick.decimal.formatter <- function(file=NULL,file_prefix,formatter="%04d",file_extension=".asc",nlayers=10,use.read.raster.from.url=FALSE,crs=NULL,start.from.zero=FALSE) {
	
	b <- NULL
	if (!is.null(file) & !is.na(file)) {
		
		file_prefix <- ""
		file_extension <- NA
		formatter <- file
		
	}
	
	
	
	startlayer <- 1  # startlayer id 1 by default but can be 0 optionally!!
	
	if (start.from.zero) startlayer <- 0
	formatter <- rep(formatter,nlayers-startlayer+1)
	for (i in startlayer:nlayers) {
		is <- i-startlayer+1
		formatter[is] <- sprintf(formatter[is],i)		
	}
	
	
	
	filepath <- paste(file_prefix,formatter,sep="")
	
	
	if (!is.na(file_extension) & !is.null(file_extension)) {
		if (str_sub(file_extension,1,1)==".")   {
			filepath <- paste(filepath,file_extension,sep="") 
		} else { 	
			filepath <- paste(filepath,file_extension,sep=".") 
		}
	}
	
	
	list <- as.list(array(NA,length(filepath)))
	names(list) <- filepath 
	
	if (use.read.raster.from.url) {
	    ####print("Warning: the option \'use.read.raster.from.url\' is not completey implemented!! Change this option if possible!")
	for (i in 1:length(list)) {
		
		
	
	
		
		list[[i]] <- read.raster.from.url(x=as.character(filepath[i]))
		
		
		}
	
	} else {
	for (i in 1:length(list)) {
		
		
		#if (!file.exists(as.character(x[i]))) print(paste("Warning Missing File:",as.character(x[i]),sep=" "))
		# read.raster.from.url
		x <- as.character(filepath[i])
		###print("list")
		####print(x)
		#####print(file.exists(x))
		if (file.exists(x)) { 
			list[[i]] <- raster(x=x)
		} else if (i>1) {
			list[[i]] <- list[[i-1]]*0.0
		}
		
		
	}
	
	}


	b <- brick(list)

	if (!is.null(crs)) {
		projection(b) <- crs
	}
	return(b)
}
	
	
