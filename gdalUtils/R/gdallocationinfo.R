#' gdallocationinfo
#' 
#' R wrapper for gdallocationinfo: raster query tool
#' 
#' @param srcfile Character. The source GDAL raster datasource name.
#' @param x Numeric. X location of target pixel. By default the coordinate system is pixel/line unless -l_srs, -wgs84 or -geoloc supplied.
#' @param y Numeric. Y location of target pixel. By default the coordinate system is pixel/line unless -l_srs, -wgs84 or -geoloc supplied.
#' @param coords Character or Matrix. Filename of coordinates (space separated, no header) or a matrix of coordinates.
#' @param xml Logical. The output report will be XML formatted for convenient post processing.
#' @param lifonly Logical. The only output is filenames production from the LocationInfo request against the database (ie. for identifying impacted file from VRT).
#' @param valonly Logical. The only output is the pixel values of the selected pixel on each of the selected bands.
#' @param b Numeric. band. Selects a band to query. Multiple bands can be listed. By default all bands are queried.
#' @param overview Numeric. overview_level. Query the (overview_level)th overview (overview_level=1 is the 1st overview), instead of the base band. Note that the x,y location (if the coordinate system is pixel/line) must still be given with respect to the base band.
#' @param l_srs Character. srs def. The coordinate system of the input x, y location.
#' @param geoloc Logical. Indicates input x,y points are in the georeferencing system of the image.
#' @param wgs84 Logical. Indicates input x,y points are WGS84 long, lat.
#' @param oo Character. "NAME=VALUE". (starting with GDAL 2.0) Dataset open option (format specific)
## @param additional_commands Character. Additional commands to pass directly to gdallocationinfo.
#' @param raw_output Logical. Dump the raw output of the gdallocationinfo (default=TRUE). If not, attempt to return a matrix of data. 
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' @return Character or matrix (if valonly=T &  raw_output=F)
#' 
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdallocationinfo' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://www.gdal.org/gdallocationinfo.html}), or,
#' in some cases, use R vectors to achieve the same end.  
#' 
#' This utility is intended to provide a variety of information about a pixel. Currently it reports three things:
#' 
#' The location of the pixel in pixel/line space.
#' The result of a LocationInfo metadata query against the datasource - currently this is only implemented for VRT files which will report the file(s) used to satisfy requests for that pixel.
#' The raster pixel value of that pixel for all or a subset of the bands.
#' The unscaled pixel value if a Scale and/or Offset apply to the band.
#' The pixel selected is requested by x/y coordinate on the commandline, or read from stdin. More than one coordinate pair can be supplied when reading coordinatesis from stdin. By default pixel/line coordinates are expected. However with use of the -geoloc, -wgs84, or -l_srs switches it is possible to specify the location in other coordinate systems.
#' 
#' The default report is in a human readable text format. It is possible to instead request xml output with the -xml switch.
#' 
#' For scripting purposes, the -valonly and -lifonly switches are provided to restrict output to the actual pixel values, or the LocationInfo files identified for the pixel.
#' 
#' It is anticipated that additional reporting capabilities will be added to gdallocationinfo in the future.
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#'
#' @references \url{http://www.gdal.org/gdallocationinfo.html}
#' 
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install
#' # and that raster and rgdal are also installed.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(valid_install)
#' {
#' 	src_dataset <- system.file("external/tahoe_highrez.tif", package="gdalUtils")
#' 	# Raw output of a single coordinate:
#' 	gdallocationinfo(srcfile=src_dataset,x=10,y=10)
#'	
#' 	# A matrix of coordinates and a clean, matrix output:
#' 	coords <- rbind(c(10,10),c(20,20),c(30,30))
#' 	gdallocationinfo(srcfile=src_dataset,coords=coords,valonly=TRUE,raw_output=FALSE)
#' }
#' @importFrom R.utils countLines
#' @importFrom utils write.table
#' @export

gdallocationinfo <- function(
		srcfile,
		x,y,
		coords,
		xml,lifonly,valonly,b,overview,l_srs,geoloc,wgs84,oo,
#		additional_commands,
		raw_output=TRUE,
		ignore.full_scan=TRUE,
		verbose=FALSE)
{	
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Start gdalinfo setup
	parameter_variables <- list(
			logical = list(
					varnames <- c("xml","lifonly","valonly",
							"geoloc","wgs84")),
			vector = list(
					varnames <- c("overview")),
			scalar = list(
					varnames <- c("x","y")),
			character = list(
					varnames <- c("l_srs","oo","srcfile")),
			repeatable = list(
					varnames <- c("b"))
	)
	
	parameter_order <- c(
			"xml","lifonly","valonly",
			"geoloc","wgs84",
			"b","overview",
			"l_srs","oo","srcfile",
			"x","y"
	)
	
	parameter_noflags <- c("srcfile","x","y")
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdallocationinfo"
	# End gdalinfo setup
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_noquotes=parameter_noquotes)
	
	# browser()
	
	if(!missing(coords))
	{
		if(class(coords)=="matrix")
		{
			# We have to write these out to disk for now, unless
			# I can figure out a more clever way of doing this...
			tempcoords_fname <- tempfile()
			write.table(coords,file=tempcoords_fname,col.names=F,row.names=F)
			cmd <- paste(cmd, "<", tempcoords_fname)
			ncoords <- nrow(coords)
		} else
		{
			# Coords is a filename:
			if(!file.exists(coords)) 
			{ 
				stop("coords must be a matrix or an existing file of coordinates.  Please fix.")
			} else
			{
				cmd <- paste(cmd, "<", coords)
				ncoords <- countLines(coords)[1]
			}
		}
	}
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	if(!missing(valonly))
	{
		if(!raw_output & valonly)
		{
			# NAs only appear once in output, need to fix:
			# browser()
			# Check for NAs:
			na_loc <- which(cmd_output == "")
			if(length(na_loc) > 0)
			{
				# We need to check for the number of bands in the input:
				if(missing(b))
				{
					nbands <- gdalinfo(srcfile,raw_output=F)$bands
				} else
				{
					nbands <- length(b)
				}
				
				missing_data_for_single_pixel <- rep(NA,nbands)
				
				# Insert function from http://stackoverflow.com/questions/18951248/insert-elements-in-a-vector-in-r
				insert.at <- function(a, pos, ...){
					dots <- list(...)
					stopifnot(length(dots)==length(pos))
					result <- vector("list",2*length(pos)+1)
					result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
					result[c(FALSE,TRUE)] <- dots
					unlist(result)
				}
				
				cmd_output <- insert.at(cmd_output,pos=na_loc,missing_data_for_single_pixel)
				cmd_output <- cmd_output[cmd_output != ""]
				
				
			}
			
			cmd_output <- matrix(as.numeric(cmd_output),nrow=ncoords,byrow=TRUE)
		}
	}
	if(verbose) { message(cmd_output) }
	
	return(cmd_output)
}
