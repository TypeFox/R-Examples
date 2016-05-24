#' gdal_grid
#' 
#' R wrapper for gdal_grid: creates regular grid from the scattered data
#' 
#' @param src_datasource Character. Any OGR supported readable datasource.
#' @param dst_filename Character. The GDAL supported output file.
#' @param ot Character. "type". For the output bands to be of the indicated data type.
#' @param of Character. "format". Select the output format. The default is GeoTIFF (GTiff). Use the short format name.
#' @param txe Numeric. c(xmin,xmax). Set georeferenced X extents of output file to be created.
#' @param tye Numeric. c(ymin,ymax). Set georeferenced Y extents of output file to be created.
#' @param outsize Numeric. c(xsize,ysize). Set the size of the output file in pixels and lines.
#' @param a_srs Character. "srs_def". Override the projection for the output file. The srs_def may be any of the usual GDAL/OGR forms, complete WKT, PROJ.4, EPSG:n or a file containing the WKT.
#' @param zfield Character. "field_name". Identifies an attribute field on the features to be used to get a Z value from. This value overrides Z value read from feature geometry record (naturally, if you have a Z value in geometry, otherwise you have no choice and should specify a field name containing Z value).
#' @param z_increase Numeric. increase_value. Addition to the attribute field on the features to be used to get a Z value from. The addition should be the same unit as Z value. The result value will be Z value + Z increase value. The default value is 0.
#' @param z_multiply Numeric. multiply_value. This is multiplication ratio for Z field. This can be used for shift from e.g. foot to meters or from elevation to deep. The result value will be (Z value + Z increase value) * Z multiply value. The default value is 1.
#' @param a Character. [algorithm[:parameter1=value1][:parameter2=value2]...] Set the interpolation algorithm or data metric name and (optionally) its parameters. See INTERPOLATION ALGORITHMS and DATA METRICS sections for further discussion of available options.
#' @param spat Numeric. c(xmin,ymin,xmax,ymax). Adds a spatial filter to select only features contained within the bounding box described by (xmin, ymin) - (xmax, ymax).
#' @param clipsrc Numeric or Character.  c(xmin,ymin,xmax,ymax)|WKT|datasource|spat_extent. Adds a spatial filter to select only features contained within the specified bounding box (expressed in source SRS), WKT geometry (POLYGON or MULTIPOLYGON), from a datasource or to the spatial extent of the -spat option if you use the spat_extent keyword. When specifying a datasource, you will generally want to use it in combination of the -clipsrclayer, -clipsrcwhere or -clipsrcsql options.
#' @param clipsrcsql Character. Select desired geometries using an SQL query instead.
#' @param clipsrclayer Character. "layername". Select the named layer from the source clip datasource.
#' @param clipsrcwhere Character. "expression". Restrict desired geometries based on attribute query.
#' @param l Character. "layername". Indicates the layer(s) from the datasource that will be used for input features. May be specified multiple times, but at least one layer name or a -sql option must be specified.
#' @param where Character. "expression". An optional SQL WHERE style query expression to be applied to select features to process from the input layer(s).
#' @param sql Character. "select_statement". An SQL statement to be evaluated against the datasource to produce a virtual layer of features to be processed.
#' @param co Character. "NAME=VALUE". Passes a creation option to the output format driver. Multiple -co options may be listed. See format specific documentation for legal creation options for each format.
#' @param q Logical. Suppress progress monitor and other non-error output.
#' @param output_Raster Logical. Return output dst_filename as a RasterBrick?
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  


#' @return NULL or if(output_Raster), a RasterBrick.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdal_grid' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdal_contour format (\url{http://www.gdal.org/gdal_grid.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' INTERPOLATION ALGORITHMS
#' 
#' There are number of interpolation algorithms to choose from.
#' \itemize{
#' \item{invdist
#' 
#' Inverse distance to a power. This is default algorithm. It has following parameters:
#' \itemize{
#' \item{power:
#' Weighting power (default 2.0).
#' }
#' \item{
#' smoothing:
#' Smoothing parameter (default 0.0).
#' }
#' \item{
#' radius1:
#' The first radius (X axis if rotation angle is 0) of search ellipse. Set this parameter to zero to use whole point array. Default is 0.0.
#' }
#' \item{
#' radius2:
#' The second radius (Y axis if rotation angle is 0) of search ellipse. Set this parameter to zero to use whole point array. Default is 0.0.
#' }
#' \item{
#' angle:
#' Angle of search ellipse rotation in degrees (counter clockwise, default 0.0).
#' }
#' \item{
#' max_points:
#' Maximum number of data points to use. Do not search for more points than this number. This is only used if search ellipse is set (both radii are non-zero). Zero means that all found points should be used. Default is 0.
#' }
#' \item{
#' min_points:
#' Minimum number of data points to use. If less amount of points found the grid node considered empty and will be filled with NODATA marker. This is only used if search ellipse is set (both radii are non-zero). Default is 0.
#' }
#' \item{
#' nodata:
#' NODATA marker to fill empty points (default 0.0).
#' }
#' }
#' }
#' \item{invdistnn
#' 
#' (Since GDAL 2.1) Inverse distance to a power with nearest neighbor searching, ideal when max_points is used. It has following parameters:
#' \itemize{
#' \item{power:
#' Weighting power (default 2.0).
#' }
#' \item{radius:
#' The radius of the search circle, which should be non-zero. Default is 1.0.
#' }
#' \item{max_points:
#' Maximum number of data points to use. Do not search for more points than this number. Found points will be ranked from nearest to furthest distance when weighting. Default is 12.
#' }
#' \item{min_points:
#' Minimum number of data points to use. If less amount of points found the grid node is considered empty and will be filled with NODATA marker. Default is 0.
#' }
#' \item{nodata:
#' NODATA marker to fill empty points (default 0.0).
#' }
#' }
#' }
#' \item{average
#' 
#' Moving average algorithm. It has following parameters:
#' \itemize{
#' \item{radius1:
#' The first radius (X axis if rotation angle is 0) of search ellipse. Set this parameter to zero to use whole point array. Default is 0.0.
#' }
#' \item{radius2:
#' The second radius (Y axis if rotation angle is 0) of search ellipse. Set this parameter to zero to use whole point array. Default is 0.0.
#' }
#' \item{angle:
#' Angle of search ellipse rotation in degrees (counter clockwise, default 0.0).
#' }
#' \item{min_points:
#' Minimum number of data points to use. If less amount of points found the grid node considered empty and will be filled with NODATA marker. Default is 0.
#' }
#' \item{nodata:
#' NODATA marker to fill empty points (default 0.0).
#' Note, that it is essential to set search ellipse for moving average method. It is a window that will be averaged when computing grid nodes values.
#' }
#' }
#' }
#' \item{nearest
#' 
#' Nearest neighbor algorithm. It has following parameters:
#' \itemize{
#' \item{radius1:
#' The first radius (X axis if rotation angle is 0) of search ellipse. Set this parameter to zero to use whole point array. Default is 0.0.
#' }
#' \item{radius2:
#' The second radius (Y axis if rotation angle is 0) of search ellipse. Set this parameter to zero to use whole point array. Default is 0.0.
#' }
#' \item{angle:
#' Angle of search ellipse rotation in degrees (counter clockwise, default 0.0).
#' }
#' \item{nodata:
#' NODATA marker to fill empty points (default 0.0).
#' }
#' }
#' }
#' \item{linear
#' 
#' (Since GDAL 2.1) Linear interpolation algorithm.
#' 
#' The Linear method performs linear interpolation by compution a Delaunay triangulation of the point cloud, finding in which triangle of the triangulation the point is, and by doing linear interpolation from its barycentric coordinates within the triangle. If the point is not in any triangle, depending on the radius, the algorithm will use the value of the nearest point or the nodata value.
#' 
#' It has following parameters:
#' \itemize{
#' \item{radius:
#' In case the point to be interpolated does not fit into a triangle of the Delaunay triangulation, use that maximum distance to search a nearest neighbour, or use nodata otherwise. If set to -1, the search distance is infinite. If set to 0, nodata value will be always used. Default is -1.
#' }
#' \item{nodata:
#' NODATA marker to fill empty points (default 0.0).
#' }
#' }
#' }
#' }
#' DATA METRICS
#' 
#' Besides the interpolation functionality gdal_grid can be used to compute some data metrics using the specified window and output grid geometry. These metrics are:
#' \itemize{
#' \item{minimum:
#' Minimum value found in grid node search ellipse.
#' }
#' \item{maximum:
#' Maximum value found in grid node search ellipse.
#' }
#' \item{range:
#' A difference between the minimum and maximum values found in grid node search ellipse.
#' }
#' \item{count:
#' A number of data points found in grid node search ellipse.
#' }
#' \item{average_distance:
#' An average distance between the grid node (center of the search ellipse) and all of the data points found in grid node search ellipse.
#' }
#' \item{average_distance_pts:
#' An average distance between the data points found in grid node search ellipse. The distance between each pair of points within ellipse is calculated and average of all distances is set as a grid node value.
#' }
#' }
#' All the metrics have the same set of options:
#' \itemize{
#' \item{radius1:
#' The first radius (X axis if rotation angle is 0) of search ellipse. Set this parameter to zero to use whole point array. Default is 0.0.
#' }
#' \item{radius2:
#' The second radius (Y axis if rotation angle is 0) of search ellipse. Set this parameter to zero to use whole point array. Default is 0.0.
#' }
#' \item{angle:
#' Angle of search ellipse rotation in degrees (counter clockwise, default 0.0).
#' }
#' \item{min_points:
#' Minimum number of data points to use. If less amount of points found the grid node considered empty and will be filled with NODATA marker. This is only used if search ellipse is set (both radii are non-zero). Default is 0.
#' }
#' \item{nodata:
#' NODATA marker to fill empty points (default 0.0).
#' }
#' }
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL that has the right drivers 
#' as specified with the "of" (output format) parameter.
#' 
#' The user can choose to (optionally) return a RasterBrick of the output file (assuming
#' raster/rgdal supports the particular output format).
#'
#' @references \url{http://www.gdal.org/gdal_grid.html}
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install
#' # and that raster and rgdal are also installed.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(require(raster) && valid_install)
#' {
#' 	# Create a properly formatted CSV:
#' 	temporary_dir <- tempdir()
#' 	tempfname_base <- file.path(temporary_dir,"dem")
#' 	tempfname_csv <- paste(tempfname_base,".csv",sep="")
#' 	
#' 	pts <- data.frame(
#' 			Easting=c(86943.4,87124.3,86962.4,87077.6),
#' 			Northing=c(891957,892075,892321,891995),
#' 			Elevation=c(139.13,135.01,182.04,135.01)
#' 	)
#' 	
#' 	write.csv(pts,file=tempfname_csv,row.names=FALSE)
#' 	
#' 	# Now make a matching VRT file
#' 	tempfname_vrt <- paste(tempfname_base,".vrt",sep="")
#' 	vrt_header <- c(
#' 	'<OGRVRTDataSource>',
#' 	'\t<OGRVRTLayer name="dem">',
#' 	'\t<SrcDataSource>dem.csv</SrcDataSource>',
#' 	'\t<GeometryType>wkbPoint</GeometryType>', 
#'  '\t<GeometryField encoding="PointFromColumns" x="Easting" y="Northing" z="Elevation"/>',
#' 	'\t</OGRVRTLayer>',
#' 	'\t</OGRVRTDataSource>'			
#' 	)
#' 	vrt_filecon <- file(tempfname_vrt,"w")
#' 	writeLines(vrt_header,con=vrt_filecon)
#' 	close(vrt_filecon)
#' 
#' 	tempfname_tif <- paste(tempfname_base,".tiff",sep="")
#' 	
#' 	# Now run gdal_grid:
#' 	setMinMax(gdal_grid(src_datasource=tempfname_vrt,
#' 		dst_filename=tempfname_tif,a="invdist:power=2.0:smoothing=1.0",
#' 		txe=c(85000,89000),tye=c(894000,890000),outsize=c(400,400),
#' 		of="GTiff",ot="Float64",l="dem",output_Raster=TRUE))
#' }
#' @import rgdal
#' @export

gdal_grid <- function(
		src_datasource,dst_filename,
		ot,of,txe,tye,outsize,a_srs,zfield,
		z_increase,z_multiply,a,spat,clipsrc,
		clipsrcsql,clipsrclayer,clipsrcwhere,
		l,where,sql,co,q,
		output_Raster=FALSE,
		ignore.full_scan=TRUE,
		verbose=FALSE #,
#		...
)
{
	if(output_Raster && (!requireNamespace("raster") || !requireNamespace("rgdal")))
	{
		warning("rgdal and/or raster not installed. Please install.packages(c('rgdal','raster')) or set output_Raster=FALSE")
		return(NULL)
	}
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Place all gdal function variables into these groupings:
	parameter_variables <- list(
			logical = list(
					varnames <- c(
							"q"
					)),
			vector = list(
					varnames <- c(
							"txe","tye","outsize","spat"
					)),
			scalar = list(
					varnames <- c(
							"z_increase","z_multiply"
					)),
			character = list(
					varnames <- c(
							"ot","of","a_srs","zfield","a","clipsrcsql","clipsrclayer",
							"l","where","sql","co","src_datasource","dst_filename"
					)),
			repeatable = list(
					varnames <- c(
					
					))
	)
	
	# Fix for clipsrc bug reported by Alex Zvoleff, 5/11/2015
	if(!missing(clipsrc))
	{
		if(is.numeric(clipsrc))
		{
			parameter_variables$vector[[1]] <- c(parameter_variables$vector[[1]],"clipsrc")
		} else
		{
			parameter_variables$character[[1]] <- c(parameter_variables$character[[1]],"clipsrc")
		}
	}
	
	parameter_order <- c(
			"q",
			"txe","tye","outsize","spat",
			"z_increase","z_multiply",
			"ot","of","a_srs","zfield","a","clipsrc","clipsrcsql","clipsrclayer",
			"l","where","sql","co","src_datasource","dst_filename"
	
	)
	
	parameter_noflags <- c("src_datasource","dst_filename")
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdal_grid"
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_noquotes=parameter_noquotes,
			#		gdal_installation_id=gdal_chooseInstallation(hasDrivers=of))
			gdal_installation_id=gdal_chooseInstallation())
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	if(output_Raster)
	{
		return(brick(dst_filename))	
	} else
	{
		return(NULL)
	}	
}