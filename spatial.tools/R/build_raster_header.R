#' Builds a raster header for a flat binary file.
#' @param x_filename Character. The filename of the input binary file.
#' @param reference_raster Raster*. A Raster* object containing the header information to be used.
#' @param out_nlayers Numeric. The number of layers in the flat binary file (defaults to nlayers(reference_raster)).
#' @param datatype Character. The dataType of the flat binary file.  See ?dataType for available datatypes.  Default is 'FLT8S'.
#' @param bandorder Character. The bandorder ('BIP','BIL','BSQ') of the file. Default is 'BSQ'.
#' @param format Character. The format of the header.  See ?hdr for valid entries.  Default is 'raster'.  CURRENTLY UNSUPPORTED.
#' @param setMinMax Logical. Set the min/max for the file (will take longer to execute)?  Default=FALSE.
#' @param additional_header Character. Create additional output headers for use with other GIS systems (see \code{\link{hdr}}). Set to NULL (default) to suppress.
#' @param verbose logical. Enable verbose execution? Default is FALSE.  
#' @author Jonathan A. Greenberg and Robert Hijimans (\email{spatial.tools@@estarcion.net})
#' @seealso \code{\link{hdr}},\code{\link{dataType}}
#' @examples \dontrun{ 
#' tahoe_highrez <- brick(system.file("external/tahoe_highrez.tif", package="spatial.tools"))
#' test_blank_file <- create_blank_raster(filename=paste(tempfile(),".gri",sep=""),
#' 	reference_raster=tahoe_highrez,nlayers=2,
#' 	create_header=FALSE,format="raster",datatype="FLT8S",bandorder="BSQ")
#' test_blank_raster <- build_raster_header(x_filename=test_blank_file,
#' 	reference_raster=tahoe_highrez,out_nlayers=2,
#' 	datatype='FLT8S',format='raster',bandorder="BSQ",setMinMax=TRUE)
#' test_blank_raster
#' }
#' @import raster
#' @export

build_raster_header <- function(x_filename,reference_raster,out_nlayers,
		datatype='FLT8S',format='raster',bandorder="BSQ",setMinMax=FALSE,
		additional_header=NULL,
		verbose=FALSE)
{
#	require("raster")
	if(missing(out_nlayers))
	{
		out_nlayers=nlayers(reference_raster)
	}
	
	if(out_nlayers==1)
	{
		outraster <- raster(reference_raster)
	} else
	{
		outraster <- brick(raster(reference_raster),nl=out_nlayers)
	}
	
	outraster@file@name <- x_filename
	outraster@file@datanotation <- datatype
	outraster@file@bandorder <- bandorder
	if(setMinMax) outraster@data@haveminmax <- TRUE	
	else outraster@data@haveminmax <- FALSE
	
	try(outhdr <- hdr(outraster, format=format),silent=TRUE)

	if(out_nlayers==1)
	{
		outraster <- raster(paste(remove_file_extension(x_filename,".gri"),".grd",sep=""))
	} else
	{
		outraster <- brick(paste(remove_file_extension(x_filename,".gri"),".grd",sep=""))
	}
	
	if(setMinMax) outraster <- setMinMax(outraster)
	else outraster@data@haveminmax <- FALSE
	
	if(!is.null(additional_header))
	{
		hdr(outraster,format=additional_header)
	}
	
	return(outraster)
}