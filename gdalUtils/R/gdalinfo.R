#' gdalinfo
#' 
#' R wrapper for gdalinfo
#' 
#' @param datasetname Character. A raster dataset name. It can be either file name.
#' @param mm Logical. Force computation of the actual min/max values for each band in the dataset?
#' @param stats Logical. Read and display image statistics. Force computation if no statistics are stored in an image.
#' @param approx_stats Logical. Read and display image statistics. Force computation if no statistics are stored in an image. However, they may be computed based on overviews or a subset of all tiles. Useful if you are in a hurry and don't want precise stats.
#' @param hist Logical. Report histogram information for all bands.
#' @param nogcp Logical. Suppress ground control points list printing. It may be useful for datasets with huge amount of GCPs, such as L1B AVHRR or HDF4 MODIS which contain thousands of them.
#' @param nomd Logical. Suppress metadata printing. Some datasets may contain a lot of metadata strings.
#' @param nrat Logical. Suppress printing of raster attribute table.
#' @param noct Logical. Suppress printing of color table.
#' @param checksum Logical. Force computation of the checksum for each band in the dataset.
#' @param mdd  Character. Report metadata for the specified domain.
#' @param nofl Logical. (GDAL >= 1.9.0) Only display the first file of the file list.
#' @param sd Numeric. (GDAL >= 1.9.0) If the input dataset contains several subdatasets read and display a subdataset with specified number (starting from 1). This is an alternative of giving the full subdataset name.
#' @param proj4 Logical. (GDAL >= 1.9.0) Report a PROJ.4 string corresponding to the file's coordinate system.
#' @param oo Character. (starting with GDAL 2.0) NAME=VALUE. Dataset open option (format specific).
#' @param version Logical. Report the version of GDAL and exit.
#' @param formats Logical. List all raster formats supported by this GDAL build (read-only and read-write) and exit. The format support is indicated as follows: 'ro' is read-only driver; 'rw' is read or write (ie. supports CreateCopy); 'rw+' is read, write and update (ie. supports Create). A 'v' is appended for formats supporting virtual IO (/vsimem, /vsigzip, /vsizip, etc). A 's' is appended for formats supporting subdatasets. Note: The valid formats for the output of gdalwarp are formats that support the Create() method (marked as rw+), not just the CreateCopy() method.
#' @param format Character. List detailed information about a single format driver. The format should be the short name reported in the --formats list, such as GTiff.
#' @param optfile Character. Read the named file and substitute the contents into the commandline options list. Lines beginning with # will be ignored. Multi-word arguments may be kept together with double quotes.
#' @param config Character. Sets the named configuration keyword to the given value, as opposed to setting them as environment variables. Some common configuration keywords are GDAL_CACHEMAX (memory used internally for caching in megabytes) and GDAL_DATA (path of the GDAL "data" directory). Individual drivers may be influenced by other configuration options.
#' @param debug Character. Control what debugging messages are emitted. A value of ON will enable all debug messages. A value of OFF will disable all debug messages. Another value will select only debug messages containing that string in the debug prefix code.
## @param additional_commands Character. Additional commands to pass directly to gdalinfo.
#' @param raw_output Logical. Dump the raw output of the gdalinfo (default=TRUE). If not, attempt to return a clean list (not all parameters will be retained, at present). 
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' 
#' @return character (if raw_output=TRUE) or list (if raw_output=FALSE).
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) and Matteo Mattiuzzi (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdalinfo' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://www.gdal.org/gdalinfo.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#' 
#' By default, this will return the gdalinfo as a character vector, one line of the output
#' per element.  The user can choose raw_output=FALSE for a cleaner format (similar to GDALinfo
#' in the rgdal package), although not all parameters are preserved.  
#'
#' @references \url{http://www.gdal.org/gdalinfo.html}
#' 
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(valid_install)
#' {
#' src_dataset <- system.file("external/tahoe_highrez.tif", package="gdalUtils")
#' # Command-line gdalinfo call:
#' # gdalinfo tahoe_highrez.tif
#' gdalinfo(src_dataset,verbose=TRUE)
#' }
#' @importFrom utils glob2rx
#' @export

gdalinfo <- function(datasetname,mm,stats,
		approx_stats,hist,nogcp,nomd,nrat,noct,nofl,checksum,
		proj4,oo,mdd,sd,
		version,formats,format,optfile,config,debug,
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
					varnames <- c("mm","stats","approx_stats","hist","nogcp","nomd",
							"nrat","noct","checksum","nofl","proj4",
							"version","formats")),
			vector = list(
					varnames <- NULL),
			scalar = list(
					varnames <- c("sd")),
			character = list(
					varnames <- c("mdd","datasetname",
							"format","optfile","config","debug","oo")),
			repeatable = list(
					varnames <- NULL)
	)
	
	parameter_order <- c(
			"mm","stats","approx_stats","hist","nogcp","nomd","nrat","noct","nofl","checksum",
			"proj4","mdd","sd",
			"version","formats","format","optfile","config","debug","oo",
			"datasetname")
	
	parameter_noflags <- c("datasetname")
	
	parameter_doubledash <- c("version","formats","format","optfile","config","debug")
	
	executable <- "gdalinfo"
	# End gdalinfo setup
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_doubledash=parameter_doubledash)
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	if(verbose) { message(cmd_output) } 
	
	# (Optional) return Raster
	if(raw_output)
	{
		return(cmd_output)	
	} else
	{
		result <- list()
		
		# browser()
		
		# Raster size (in pixels and lines).
		dims          <- strsplit(gsub(grep(cmd_output,pattern="Size is ",value=TRUE), pattern="Size is ",replacement=""),",")[[1]]
		result$rows   <- as.numeric(dims[2])
		result$columns<- as.numeric(dims[1])
		
		# Bands!
		bands <- grep(cmd_output,pattern="Band ",value=TRUE)
		if(length(bands)==0) result$bands=1 else
		{
			result$bands <- length(bands)
		}
		
		orig          <- as.numeric(strsplit(gsub(strsplit(grep(cmd_output,pattern="Lower Left  \\(",value=TRUE), "\\) \\(")[[1]][1],pattern="Lower Left  \\(",replacement=""),",")[[1]])
		result$ll.x   <- orig[1]
		result$ll.y   <- orig[2]
		
		res           <- as.numeric(strsplit(gsub(gsub(grep(cmd_output,pattern="Pixel Size = \\(",value=TRUE),pattern="Pixel Size = \\(",replacement=""),pattern="\\)",replacement=""),",")[[1]])
		result$res.x  <- res[1]
		result$res.y  <- res[2]
		
		result$file <- gsub(grep(cmd_output,pattern="Files: ",value=TRUE),pattern="Files: ",replacement="")
		
		if(!missing(proj4))
		{
			# Thanks to http://haotu.wordpress.com/2011/12/02/trim-remove-trailing-and-leading-spaces-from-a-character-string-in-r/
			# for hints.
			if(proj4) result$proj4 <- sub("\\s+$","",gsub("'","",cmd_output[grep(pattern="PROJ.4 string is:",cmd_output)+1]))			
		}
#		result$oblique.x <- NA
#		result$oblique.y <- NA
		
		# The format driver used to access the file.
		result$driver <- strsplit(gsub(grep(cmd_output,pattern="Driver: ",value=TRUE), pattern="Driver: ",replacement=""),"/")[[1]][1]
		
		# The coordinate system for the file (in OGC WKT).
		# TODO
		
		# The geotransform associated with the file (rotational coefficients are currently not reported).
		# TODO
		
		# Corner coordinates in georeferenced, and if possible lat/long based on the full geotransform (but not GCPs).
		ul <- as.numeric(strsplit(strsplit(strsplit(cmd_output[grep(pattern=glob2rx("Upper Left*"),cmd_output)],"\\(")[[1]][2],"\\)")[[1]][1],",")[[1]])
		ll <- as.numeric(strsplit(strsplit(strsplit(cmd_output[grep(pattern=glob2rx("Lower Left*"),cmd_output)],"\\(")[[1]][2],"\\)")[[1]][1],",")[[1]])
		ur <- as.numeric(strsplit(strsplit(strsplit(cmd_output[grep(pattern=glob2rx("Upper Right*"),cmd_output)],"\\(")[[1]][2],"\\)")[[1]][1],",")[[1]])
		lr <- as.numeric(strsplit(strsplit(strsplit(cmd_output[grep(pattern=glob2rx("Lower Right*"),cmd_output)],"\\(")[[1]][2],"\\)")[[1]][1],",")[[1]])
		corners_rbind <- rbind(ul,ll,ur,lr)
		
		result$bbox <- matrix(c(min(corners_rbind[,1]),max(corners_rbind[,1]),min(corners_rbind[,2]),max(corners_rbind[,2])),nrow=2,ncol=2,byrow=TRUE)
		colnames(result$bbox) <- c("min","max")
		rownames(result$bbox) <- c("s1","s2")
	
		# Ground control points.
		# TODO
		
		# File wide (including subdatasets) metadata.
		# TODO
		
		# Band data types.
		# TODO
		
		# Band color interpretations.
		# TODO	
		
		# Band block size.
		# TODO
		
		# Band descriptions.
		# TODO
		
		# Band min/max values (internally known and possibly computed).
		# TODO
		
		# Band checksum (if computation asked).
		# TODO
		
		# Band NODATA value.
		# TODO
		
		# Band overview resolutions available.
		# TODO
		
		# Band unit type (i.e.. "meters" or "feet" for elevation bands).
		# TODO
		
		# Band pseudo-color tables.
		# TODO	
		
		return(result)
		
	}	
}
