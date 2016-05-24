#' ogrtindex
#' 
#' R wrapper for ogrtindex: creates a tileindex
#' 
#' @param output_dataset Character. Output tile index.
#' @param src_dataset Character. Input geospatial files.
#' @param lnum Numeric. n. Add layer number 'n' from each source file in the tile index.
#' @param lname Character. name. Add the layer named 'name' from each source file in the tile index.
#' @param f Character. output_format. Select an output format name. The default is to create a shapefile.
#' @param tileindex Character. file_name. The name to use for the dataset name. Defaults to LOCATION.
#' @param write_absolute_path Logical. Filenames are written with absolute paths.
#' @param skip_different_projection Logical. Only layers with same projection ref as layers already inserted in the tileindex will be inserted.
#' @param accept_different_schemas Logical. By default ogrtindex checks that all layers inserted into the index have the same attribute schemas. If you specify this option, this test will be disabled. Be aware that resulting index may be incompatible with MapServer!
#' @param output_Vector Logical. Return output output_dataset as a Spatial* object.  Currently only works with f="ESRI Shapefile".
## @param additional_commands Character. Additional commands to pass directly to ogrtindex.
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
## @param ... Other parameters to pass to gdal_translate.


#' @return NULL or SpatialPolygonsDataFrame
#' 
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'ogrtindex' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://gdal.org/ogrtindex.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' The ogrtindex program can be used to create a tileindex - a file containing a 
#' list of the identities of a bunch of other files along with there spatial 
#' extents. This is primarily intended to be used with MapServer for tiled access 
#' to layers using the OGR connection type.
#' 
#' If no -lnum or -lname arguments are given it is assumed that all layers in 
#' source datasets should be added to the tile index as independent records.
#' 
#' If the tile index already exists it will be appended to, otherwise it will be 
#' created.
#' 
#' It is a flaw of the current ogrtindex program that no attempt is made to copy 
#' the coordinate system definition from the source datasets to the tile index 
#' (as is expected by MapServer when PROJECTION AUTO is in use).
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#'
#' @references \url{http://www.gdal.org/ogrtindex.html}
#' 
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(require(rgdal) && valid_install)
#' {
#' 	tempindex <- tempfile(fileext=".shp")
#' 	src_dir <- system.file("external/", package="gdalUtils")
#' 	src_files <- list.files(src_dir,pattern=".shp",full.names=TRUE)
#' 	ogrtindex(output_dataset=tempindex,src_dataset=src_files,
#' 			accept_different_schemas=TRUE,output_Vector=TRUE)
#' }
#' @import rgdal
#' @export

ogrtindex <- function(
		output_dataset,src_dataset,
		lnum,lname,f,tileindex,write_absolute_path,
		skip_different_projection,accept_different_schemas,		
#		additional_commands,
		output_Vector=FALSE,
		ignore.full_scan=TRUE,
		verbose=FALSE#,
#		...
		)
{
	
	if(output_Vector && !requireNamespace("rgdal"))
	{
		warning("rgdal not installed. Please install.packages('rgdal') or set output_Vector=FALSE")
		return(NULL)
	}
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Start gdalinfo setup
	parameter_variables <- list(
			logical = list(
					varnames <- c("write_absolute_path",
							"skip_different_projection",
							"accept_different_schemas")),
			vector = list(
					varnames <- c("lnum")),
			scalar = list(
					varnames <- c()),
			character = list(
					varnames <- c("output_dataset",
							"lname","f","tileindex")),
			repeatable = list(
					varnames <- c("src_dataset"))
	)
	
	# browser()
	
	parameter_order <- c(
			"write_absolute_path",
			"skip_different_projection",
			"accept_different_schemas",
			"lnum",
			"lname","f","tileindex",
			"output_dataset","src_dataset"
			)
	
	parameter_noflags <- c("output_dataset","src_dataset")
	
	parameter_doubledash <- NULL
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "ogrtindex"
	# End gdalinfo setup
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_doubledash=parameter_doubledash,
			parameter_noquotes=parameter_noquotes)
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	if(output_Vector)
	{
		return(
				readOGR(dsn=dirname(output_dataset),
						layer=basename(remove_file_extension(output_dataset))))
	} else
	{
		return(NULL)
	}
}
