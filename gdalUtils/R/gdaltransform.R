#' gdaltransform
#' 
#' R wrapper for gdaltransform: transforms coordinates
#' 
#' @param srcfile Character. File with source projection definition or GCP's. If not given, source projection is read from the command-line -s_srs or -gcp parameters.
#' @param dstfile Character. File with destination projection definition.
#' @param coords Matrix. A two-column matrix with coordinates.
#' @param s_srs Character. source spatial reference set. The coordinate systems that can be passed are anything supported by the OGRSpatialReference.SetFromUserInput() call, which includes EPSG PCS and GCSes (ie. EPSG:4296), PROJ.4 declarations (as above), or the name of a .prf file containing well known text.
#' @param t_srs Character. target spatial reference set. The coordinate systems that can be passed are anything supported by the OGRSpatialReference.SetFromUserInput() call, which includes EPSG PCS and GCSes (ie. EPSG:4296), PROJ.4 declarations (as above), or the name of a .prf file containing well known text.
#' @param to Character. "NAME=VALUE". set a transformer option suitable to pass to GDALCreateGenImgProjTransformer2(). 
#' @param order Numeric. order of polynomial used for warping (1 to 3). The default is to select a polynomial order based on the number of GCPs.
#' @param tps Logical. Force use of thin plate spline transformer based on available GCPs.
#' @param rpc Logical. Force use of RPCs.
#' @param geoloc Logical. Force use of Geolocation Arrays.
#' @param i Logical. Inverse transformation: from destination to source.
#' @param gcp Character. pixel line easting northing [elevation]: Provide a GCP to be used for transformation (generally three or more are required)
#' @param output_xy Logical. (GDAL >= 2.0) Restrict output to "x y" instead of "x y z"
## @param additional_commands Character. Additional commands to pass directly to gdaladdo.
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
## @param ... Other parameters to pass to gdaltransform.

#' @return Numeric.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdaltransform' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://gdal.org/gdaltransform.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#'
#' @references \url{http://www.gdal.org/gdaltransform.html}
#' 
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(valid_install)
#' {
#' pts <- matrix(c(177502,311865,177503,311866),ncol=2,byrow=TRUE)
#' gdaltransform(s_srs="EPSG:28992",t_srs="EPSG:31370",coords=pts,verbose=TRUE)
#' }
#' @importFrom utils write.table
#' @export

gdaltransform <- function(srcfile,dstfile,
		coords,
		s_srs,t_srs,to,order,tps,rpc,geoloc,i,gcp,output_xy,
#		additional_commands,
		ignore.full_scan=TRUE,
		verbose=FALSE#,
#		...
)
{
	
	if(missing(coords)) stop("You must set coords...")
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Start gdalinfo setup
	parameter_variables <- list(
			logical = list(
					varnames <- c("tps","rpc","geoloc","i","output_xy"
					)),
			vector = list(
					varnames <- c()),
			scalar = list(
					varnames <- c("order")),
			character = list(
					varnames <- c("s_srs","t_srs","to",
							"srcfile","dstfile")),
			repeatable = list(
					varnames <- c("gcp"))
	)
	
	parameter_order <- c(
			"tps","rpc","geoloc","i","output_xy","order","s_srs","t_srs","to",
			"srcfile","dstfile")
	
	parameter_noflags <- c("srcfile","dstfile")
	
	parameter_doubledash <- NULL
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdaltransform"
	# End gdalinfo setup
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_doubledash=parameter_doubledash,
			parameter_noquotes=parameter_noquotes)
	
	# gdaltransform defaults to reading from the command line, so we'll
	# hack this and make a temp file for it to read.
	tempcoords_fname <- tempfile()
	write.table(coords,file=tempcoords_fname,col.names=F,row.names=F)
	cmd <- paste(cmd, "<", tempcoords_fname)
	
	# browser()
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	# Convert to numeric:
	
	output_coords <- matrix(as.numeric(unlist(strsplit(cmd_output,split=" "))),ncol=3,byrow=T)
	
	return(output_coords)
}
