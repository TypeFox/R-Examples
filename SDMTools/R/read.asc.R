#' ESRI ASCII Raster File Import And Export
#' 
#' \code{read.asc} and \code{read.asc.gz} reads ESRI ArcInfo ASCII raster file
#' either uncompressed or compressed using gzip. \cr \cr \code{write.asc} and
#' \code{write.asc.gz} writes an asc object to a ESRI ArcInfo ASCII raster
#' file. The output can be either compressed or uncompressed. \cr \cr These
#' functions are faster methods based on the adehabitat import.asc and
#' export.asc.\cr\cr \code{write.asc2} and \code{write.asc2.gz} are even faster
#' implementations but have less error checking. \cr \cr \code{image.asc} and
#' \code{print.asc} are generic methods associated with plotting & summarizing
#' data of class 'asc'; they were modified from adehabitat package.
#' 
#' Implements a faster version of import.asc or export.asc from the adehabitat
#' package. In addition, files can be read in and written to in gzip compressed
#' format.\cr\cr Generic methods of print and image were modified from
#' adehabitat. Further details of them are found there.
#' 
#' @param file a character string representing the filename of the input/output
#' file. The file extension should always be '.asc'.
#' @param gz defines if the object is or should be compressed using gzip
#' @param x an object of class 'asc' as defined in the adehabitat package
#' @param sigdig is the number of significant digits to write when creating the
#' ascii grid file
#' @param col for maps of type \code{"numeric"}, the colors to be used (see
#' \code{help(par)})
#' @param clfac for maps of type \code{"factor"}, a character vector giving the
#' names of colors for each level of the factor (see \code{help(colasc)})
#' @param \dots additional arguments to be passed to the generic function
#' \code{image} or \code{print}
#' @return Returns a raster matrix of the class 'asc' defined in the adehabitat
#' package with the following attributes: \item{xll}{the x coordinate of the
#' center of the lower left pixel of the map} \item{yll}{the y coordinate of
#' the center of the lower left pixel of the map} \item{cellsize}{the size of a
#' pixel on the studied map} \item{type}{either 'numeric' or 'factor'}
#' \item{levels}{if type = 'factor', the levels of the factor.}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' 
#' 
#' #create a simple object of class 'asc'
#' tasc = as.asc(matrix(rep(x=1:10, times=1000),nr=100)); print(tasc)
#' 
#' #write out the raster grid file
#' write.asc(tasc,'t.raster.asc')
#' write.asc.gz(tasc,'t.raster.asc') #actually save file name as t.raster.asc.gz
#' 
#' #read in the raster grid files
#' tasc2 = read.asc('t.raster.asc')
#' tasc3 = read.asc.gz('t.raster.asc.gz')
#' 
#' #remove the temporary raster
#' unlink(c('t.raster.asc','t.raster.asc.gz'))
#' 
#' 
#' @export 
#' @useDynLib SDMTools writeascdata
read.asc <-
function (file, gz=FALSE) {
  #confirm ascii grid file name is specified
	if (gz) {
		if (substr(file, nchar(file) - 2, nchar(file)) != ".gz") stop("not a valid .gz file")
	} else {
		if (substr(tolower(file), nchar(file) - 3, nchar(file)) != ".asc") stop("not a valid .asc file")
	}
	
	#read in the header
	if (gz) { zz <- gzfile(file, "r") } else { zz <- file(file, "r") }
		nc <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); nc <- as.numeric(nc[[2]][1])#number of columns
		nl <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); nl <- as.numeric(nl[[2]][1])#number of rows
		xll <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); #lower left corner
		yll <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); #lower left corner
		cs <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); cs <- as.numeric(cs[[2]][1])#cell size
		nas <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); nas <- as.numeric(nas[[2]][1])#nodata value
	#close the link to the file
	close(zz)
	
	#ensure xll & yll are centers of the cells
	if ((xll[[1]][1] == "xllcenter") | (xll[[1]][1] == "XLLCENTER")) { xll=as.numeric(xll[[2]][1]) } else { xll=as.numeric(xll[[2]][1])+ cs/2 }
	if ((yll[[1]][1] == "yllcenter") | (xll[[1]][1] == "YLLCENTER")) { yll=as.numeric(yll[[2]][1]) } else { yll=as.numeric(yll[[2]][1])+ cs/2 }

	#read in the data skipping the first six header rows
	if (gz) { zz <- gzfile(file, "r"); output <- scan(zz,nmax=nl*nc,skip=6,quiet = TRUE); close(zz);
	} else { output <- scan(file,nmax=nl*nc,skip=6, quiet = TRUE) }

	#convert no data to NA
	output[output == nas] <- NA
	#convert data to matrix
	output <- matrix(c(as.matrix(output)), ncol = nl, nrow = nc)
	output <- output[, ncol(output):1]
	#define the attributes
	attr(output, "xll") <- xll
	attr(output, "yll") <- yll
	attr(output, "cellsize") <- cs
	attr(output, "type") <- 'numeric'
	class(output) <- "asc"
	#return the file
	return(output)
}

#' @rdname read.asc
#' @export
read.asc.gz <-
function (file) {
	return(read.asc(file, gz=TRUE))
}

#' @rdname read.asc
#' @export
write.asc <-
function (x, file, gz=FALSE) {
  #confirm asc object and file named appropriately
  if (!inherits(x, "asc")) stop("Non convenient data")
  if (substr(file, nchar(file) - 3, nchar(file)) != ".asc") file <- paste(file, ".asc", sep = "")
  #open a connection to file
  if (gz) { zz <- gzfile(paste(file, ".gz", sep = ""),"w") } else { zz <- file(file, "w") }
    #write the header info
    cat("ncols         ",nrow(x),'\n',sep = "",file=zz)
    cat("nrows         ",ncol(x),'\n',sep = "",file=zz)
    cat("xllcorner     ",as.character(attr(x,"xll")-attr(x,"cellsize")/2),'\n',sep = "",file=zz)
    cat("yllcorner     ",as.character(attr(x,"yll")-attr(x,"cellsize")/2),'\n',sep = "",file=zz)
    cat("cellsize      ",as.character(attr(x, "cellsize")),'\n',sep = "",file=zz)
    cat("NODATA_value  ", -9999,'\n',sep = "",file=zz)
    #prep and write the data
    x[is.na(x)] <- -9999 #change na values
    x <- x[, ncol(x):1] #reorder
    x <- do.call('rbind',list(x,"\n")) #add new line character
    cat(x,file=zz)
  #close the connection to the file
  close(zz)
}

#' @rdname read.asc
#' @export
write.asc.gz <-
function (x, file) {
	write.asc(x, file, gz=TRUE)
}

#' @rdname read.asc
#' @export
#' @import R.utils
write.asc2 <-
function (x, file, sigdig = 0, gz=FALSE) {
	###confirm asc object and file named appropriately
	if (!inherits(x, "asc")) stop("Non convenient data")
	if (substr(file, nchar(file) - 3, nchar(file)) != ".asc") file <- paste(file, ".asc", sep = "")
	###write out the data
	tt = .Call('writeascdata', ncol(x) , nrow(x) , as.character(attr(x,"xll")-attr(x,"cellsize")/2) ,
		as.character(attr(x,"yll")-attr(x,"cellsize")/2) , as.character(attr(x, "cellsize")) , x , file , sigdig, PACKAGE='SDMTools')
	if (gz) {
		gzip(file)
	}
}

#' @rdname read.asc
#' @export
#' @import R.utils
write.asc2.gz <-
function (x, file, sigdig = 0) {
	write.asc2(x, file, sigdig=sigdig, gz=TRUE)
}

#' @rdname read.asc
#' @export
"image.asc" <- function (x, col = gray((240:1)/256), clfac = NULL, ...)
{
    ## Verifications
    if (!inherits(x, "asc")) stop("not an \"asc\" object")

    ## Coordinates of the pixels
    xy <- getXYcoords(x)
    xx <- xy$x
    yy <- xy$y

    ## If the variable is numeric
    if (attr(x, "type") == "numeric") image(x = xx, y = yy, x, asp = 1, col = col, ...)

    ## For a factor: creates colors
    if (attr(x, "type") == "factor") {
        if (is.null(clfac)) {
            clfac <- rainbow(nlevels(x))
            clfac <- clfac[as.numeric(levels(factor(x)))]
        }
        image(x = xx, y = yy, x, asp = 1, col = clfac, ...)
    }
}

#' @rdname read.asc
#' @export
"print.asc" <- function(x, ...)
{
    ## Verifications
    if (!inherits(x, "asc")) stop("Non convenient data")

    ## The output
    cat("Raster map of class \"asc\":\n")
    cat("Cell size: ", attr(x, "cellsize"), "\n")
    cat("Number of rows: ", ncol(x), "\n")
    cat("Number of columns: ", nrow(x), "\n")
    cat("Type: ", attr(x, "type"), "\n")
}

