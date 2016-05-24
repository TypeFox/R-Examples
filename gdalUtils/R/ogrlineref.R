#' ogrlineref
#' 
#' R wrapper for ogrlineref: create a linear reference
#' 
#' @param help_general Logical. Show the usage.
#' @param progress Logical. Show progress.
#' @param quiet Logical. Suppress all messages except errors and results.
#' @param f Character. format_name. Select an output format name. The default is to create a shapefile.
#' @param dsco Character. "NAME=VALUE". Dataset creation option (format specific).
#' @param lco Character. "NAME=VALUE". Layer creation option (format specific).
#' @param create Logical. Create the linear reference file (linestring of parts).
#' @param l Character. src_line_datasource_name. The path to input linestring datasource (e.g. the road)
#' @param ln Character. layer_name. The layer name in datasource
#' @param lf Character. field_name. The field name of uniq values to separate the input lines (e.g. the set of roads)
#' @param p Character. src_repers_datasource_name. The path to linear references points (e.g. the road mile-stones)
#' @param pn Character. layer_name. The layer name in datasource.
#' @param pm Character. pos_field_name.The field name of distances along path (e.g. mile-stones values)
#' @param pf Character. field_name. The field name of uniq values to map input reference points to lines
#' @param r Character. src_parts_datasource_name. The path to linear reference file
#' @param rn Character. layer_name. The layer name in datasource.
#' @param o Character. dst_datasource_name. The path to output linear reference file (linestring datasource)
#' @param on Character. layer_name. The layer name in datasource.
#' @param of Character. field_name. The field name for storing the uniq values of input lines
#' @param s Numeric. step. The part size in linear units.
#' @param get_pos Logical. Return linear referenced postion for input X, Y
#' @param x Numeric. long. Input X coordinate
#' @param y Numeric. lat. Input Y coordinate
#' @param get_coord Logical. Return point on path for input linear distance.
#' @param m Numeric. position. The input linear distance
#' @param get_subline Logical. Return the portion of the input path from and to input linear positions
#' @param mb Numeric. position. The input begin linear distance.
#' @param me Numeric. position. The input end linear distance
## @param additional_commands Character. Additional commands to pass directly to ogrlineref.
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
## @param ... Other parameters to pass to ogrlineref.
#' 
#' @return NULL
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'ogrlineref' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://gdal.org/ogrlineref.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' The utility can be used for:
#' \itemize{
#' \item{create linear reference file from input data}
#' \item{return the "linear referenced" distance for the projection of the input coordinates (point) on the path}
#' \item{return the coordinates (point) on the path according to the "linear referenced" distance}
#' \item{return the portion of the path according to the "linear referenced" begin and end distances}
#' }
#' The ogrlineref program can be used to create a linear reference - a file containing 
#' a segments of special length (e.g. 1 km in reference units) and get coordinates, 
#' linear referenced distances or sublines (subpaths) from this file. The utility not 
#' required the M or Z values in geometry. The results can be stored in any OGR supported
#' format. Also some information writed to the stdout.
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#'
#' @references \url{http://www.gdal.org/ogrlineref.html}
#' 
#' @examples 
#' # No examples ATM for this function.
#' @export

ogrlineref <- function(
		help_general,
		progress,quiet,	
		f,dsco,lco,create,
		l,ln,lf,p,pn,pm,pf,
		r,rn,o,on,of,s,
		get_pos,x,y,get_coord,m,
		get_subline,mb,me,		
#		additional_commands,
		ignore.full_scan=TRUE,
		verbose=FALSE#,
#		...
		)
{
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Start gdalinfo setup
	parameter_variables <- list(
			logical = list(
					varnames <- c("help_general",
							"progress","quiet","create",
							"get_pos",
							"get_coord",
							"get_subline")),
			vector = list(
					varnames <- c("s","x","y","m","mb","me")),
			scalar = list(
					varnames <- c()),
			character = list(
					varnames <- c("f","dsco","ls","ln","lf","p","pn",
							"pm","pf","r","rn","o","on","of")),
			repeatable = list(
					varnames <- c())
	)
	
	parameter_order <- c(
			"help_general",
			"progress","quiet","create",
			"get_pos",
			"get_coord",
			"get_subline",
			"s","x","y","m","mb","me",
			"f","dsco","ls","ln","lf","p","pn",
			"pm","pf","r","rn","o","on","of"
			)
	
	parameter_noflags <- c()
	
	parameter_doubledash <- NULL
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "ogrlineref"
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
	
	return(NULL)
}
