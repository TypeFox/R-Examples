#' gdaltindex
#' 
#' R wrapper for gdaltindex: Builds a shapefile as a raster tileindex
#' 
#' @param index_file Character. The name of the output file to create/append to. The default shapefile will be created if it doesn't already exist, otherwise it will append to the existing file.
#' @param gdal_file Character. The input GDAL raster files, can be multiple files separated by spaces. Wildcards my also be used. Stores the file locations in the same style as specified here, unless -write_absolute_path option is also used.
#' @param f Character. format. The OGR format of the output tile index file. Default is Esri Shapefile.
#' @param tileindex Character. field_name. The output field name to hold the file path/location to the indexed rasters. The default tile index field name is location.
#' @param write_absolute_path Logical. The absolute path to the raster files is stored in the tile index file. By default the raster filenames will be put in the file exactly as they are specified on the command line.
#' @param skip_different_projection Logical. Only files with same projection as files already inserted in the tileindex will be inserted (unless -t_srs is specified). Default does not check projection and accepts all inputs.
#' @param t_srs Character. target_srs. Geometries of input files will be transformed to the desired target coordinate reference system. Using this option generates files that are not compatible with MapServer < 6.4. Default creates simple rectangular polygons in the same coordinate reference system as the input rasters.
#' @param src_srs_name Character. field_name. The name of the field to store the SRS of each tile. This field name can be used as the value of the TILESRS keyword in MapServer >= 6.4.
#' @param src_srs_format Character. type. The format in which the SRS of each tile must be written. Types can be AUTO, WKT, EPSG, PROJ.
#' @param lyr_name Character. name. Layer name to create/append to in the output tile index file.
#' @param output_Vector Logical. Return output dst_filename as a Spatial* object.  Currently only works with f="ESRI Shapefile".
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  

#' @return NULL or if(output_Vector), a SpatialPolygonsDataFrame.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdaltindex' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdaltindex format (\url{http://www.gdal.org/gdaltindex.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL that has the right drivers 
#' as specified with the "of" (output format) parameter.
#'  
#' The user can choose to (optionally) return a SpatialPolygonsDataFrame of the output file.
#'
#' @references \url{http://www.gdal.org/gdaltindex.html}
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install
#' # and that raster and rgdal are also installed.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(require(rgdal) && valid_install)
#' {
#' # Modified example from the original gdaltindex documentation:
#' src_folder <- system.file("external/", package="gdalUtils")
#' output_shapefile <- paste(tempfile(),".shp",sep="")
#' # Command-line gdalwarp call:
#' # gdaltindex doq_index.shp external/*.tif
#' gdaltindex(output_shapefile,list.files(path=src_folder,pattern=glob2rx("*.tif"),full.names=TRUE),
#' 	output_Vector=TRUE,verbose=TRUE)
#' }
#' @import rgdal
#' @export

gdaltindex <- function(
		index_file,gdal_file,
		f,tileindex,write_absolute_path,
		skip_different_projection,
		t_srs,src_srs_name,
		src_srs_format,
		lyr_name,
		output_Vector=FALSE,
		ignore.full_scan=TRUE,
		verbose=FALSE)
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
	
	# Place all gdal function variables into these groupings:
	parameter_variables <- list(
			logical = list(
					varnames <- c(
							"write_absolute_path",
							"skip_different_projection"
					)),
			vector = list(
					varnames <- c(
					)),
			scalar = list(
					varnames <- c(
							
					)),
			character = list(
					varnames <- c(
							"f","tileindex","t_srs",
							"src_srs_name","src_srs_format",
							"lyr_name","index_file",
							"gdal_file"
							
					)),
			repeatable = list(
					varnames <- c(
					))
	)
	
	parameter_order <- c(
			"write_absolute_path",
			"skip_different_projection",
			"f","tileindex","t_srs",
			"src_srs_name","src_srs_format",
			"lyr_name","index_file",
			"gdal_file"
	)
	
	parameter_noflags <- c("index_file","gdal_file")
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdaltindex"
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_noquotes=parameter_noquotes,
			gdal_installation_id=gdal_chooseInstallation())
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	if(output_Vector)
	{
		return(
				readOGR(dsn=dirname(index_file),
						layer=basename(remove_file_extension(index_file))))
	} else
	{
		return(NULL)
	}		
}