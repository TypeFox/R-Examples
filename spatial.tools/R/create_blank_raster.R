#' Create an empty raster and header.
#' 
#' This function creates an arbitrarily large raster
#' as a flat binary file with (optionally) a proper header 
#' for use with other functions.  This should create the blank 
#' files very quickly, as it is using some OS tricks
#' to carve out a block of space rather than writing
#' a bunch of 0s to disk sequentially.
#' 
#' @param filename Character. The output base filename of the blank file.  Will use tempfile() if nothing is provided.
#' @param format Character.  Output format.  Currently only supports "raster".
#' @param datatype Character.  Output number type.  See ?dataType.  Default is "FLT8S".  
#' @param bandorder Character.  Output band interleave.  Currently only supports "BSQ".
#' @param nrow Numeric. Number of rows of the output raster. Defaults to nrow(reference_raster).
#' @param ncol Numeric. Number of columns of the output raster. Defaults to ncol(reference_raster).
#' @param nlayers Numeric. Number of layers of the output raster Defaults to nlayer(reference_raster).
#' @param create_header Logical. Create a properly formatted header for the blank file?
#' @param reference_raster Raster*. Reference raster to derive other information, e.g. resolution, projection, datum.
#' @param return_filename Logical. Return filename of the binary file (if TRUE, default) or the Raster* itself (if FALSE).
#' @param additional_header Character. Create additional output headers for use with other GIS systems (see \code{\link{hdr}}). Set to NULL (default) to suppress.
#' @param overwrite Logical. Overwrite an existing file with the same name?
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' @return A character vector (return_filename==TRUE) or as Raster* object (return_filename==TRUE)
#' @seealso \code{\link{hdr}}
#' @details create_blank_raster is designed to quickly create a binary file of the appropriate
#' size using tricks with seek()/writeBin(). A large file can be created in a fraction of a second
#' using this function. This file could, for example, be used with mmap to realize asynchronous
#' or, OS permitting, parallel writes to a single file.  Note that setMinMax are NOT performed
#' on the output file (to save time), so on some systems you may see a warning.
#' 
#' Binary files of this type are used by a number of raster formats, including raster and ENVI.
#' 
#' @author Jonathan A. Greenberg
#' @examples \dontrun{ 
#' tahoe_highrez <- brick(system.file("external/tahoe_highrez.tif", package="spatial.tools"))
#' test_blank_file <- create_blank_raster(reference_raster=tahoe_highrez)
#' file.info(test_blank_file)
#' test_blank_raster <- create_blank_raster(reference_raster=tahoe_highrez,return_filename=FALSE)
#' test_blank_raster
#' }
#' @import raster
#' @export

create_blank_raster <- function(filename=NULL,
	format="raster",datatype="FLT8S",bandorder="BSQ",
	nrow=NULL,ncol=NULL,nlayers=NULL,
	create_header=TRUE,reference_raster=NULL,
	return_filename=TRUE,
	additional_header=NULL,
	overwrite=FALSE,verbose=FALSE)
{
	if(is.null(nrow)&!is.null(reference_raster)) nrow=nrow(reference_raster)
	if(is.null(ncol)&!is.null(reference_raster)) ncol=ncol(reference_raster)
	if(is.null(nlayers)&!is.null(reference_raster)) nlayers=nlayers(reference_raster)
	
	# Setup blank file.
	outdata_ncells=as.numeric(nrow)*as.numeric(ncol)*as.numeric(nlayers)
	if(verbose) cat("outdata_ncells=",outdata_ncells,"\n")
	if(is.null(filename))
	{	
		filename <- tempfile()
		if(verbose) { message(paste("No output file given, using a tempfile name:",filename,sep=" ")) }
		if(!file.exists(tempdir())) dir.create(tempdir())
	} 

	numBytes = dataSize(datatype)
	
#	numBytes = substr(dataType,4,4)
	
#	if(dataType=="FLT4S") numBytes = 4 # Single Precision Floating Point
#	if(dataType=="FLT8S") numBytes = 8 # Double Precision Floating Point
	
	# I've been warned about using seek on Windows, but this appears to work...
	if(verbose) { message("Creating empty file.") }
	out=file(filename,"wb")
	seek(out,(outdata_ncells-1)*numBytes)
	writeBin(raw(numBytes),out)
	close(out)
	
	if(create_header)
	{
		# Setup header.
		if(verbose) { message("Setting up output header.") }
		if(nlayers > 1) 
		{ 
#			if(nlayers(reference_raster) > 1) { 
			reference_raster=brick(raster(reference_raster),nl=nlayers) 
#			} else
#			{
#				reference_raster <- brick(raster(reference_raster),nl=nlayers) 
#			}
		} else
		{
			if(nlayers(reference_raster) > 1) { 
				reference_raster=raster(reference_raster,layer=1) 
			} 	
		}
		
		if(format=="raster") { 
			if(verbose) { message("Outformat is raster.  Appending .gri to filename.") }
			file.rename(filename,paste(filename,".gri",sep=""))
			filename=paste(filename,".gri",sep="")
		}
		
		outraster <- build_raster_header(x_filename=filename,
				reference_raster=reference_raster,
				out_nlayers=nlayers,datatype=datatype,format=format,
				bandorder=bandorder,additional_header=additional_header)
	}
	if(return_filename) return(filename) else
		return(outraster)
}