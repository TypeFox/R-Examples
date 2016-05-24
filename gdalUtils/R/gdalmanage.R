#' gdalmanage
#' 
#' R wrapper for gdalmanage: Identify, delete, rename and copy raster data files
#' 
#' @param mode Character. Mode of operation. "identify" | "copy" | "rename" | "delete". See details.
#' @param r Logical. Recursively scan files/folders for raster files.
#' @param u Logical. Report failures if file type is unidentified.
#' @param f Character. format. Specify format of raster file if unknown by the application. Uses short data format name (e.g. GTiff).
#' @param datasetname Character. Raster file to operate on.
#' @param newdatasetname Character. For copy and rename modes, you provide a source filename and a target filename, just like copy and move commands in an operating system.
## @param additional_commands Character. Additional commands to pass directly to gdaladdo.
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
## @param ... Other parameters to pass to gdaladdo.
#' 
#' @return Character.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdalmanage' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://gdal.org/gdalmanage.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' Mode of operation
#' \itemize{
#' \item{mode="identify",datasetname: List data format of file.}
#' \item{mode="copy",datasetname,newdatasetname: Create a copy of the raster file with a new name.}
#' \item{mode="rename",datasetname,newdatasetname: Change the name of the raster file.}
#' \item{mode="delete",datasetname: Delete raster file.}
#' }
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#'
#' @references \url{http://www.gdal.org/gdalmanage.html}
#' 
#' @examples 
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(valid_install)
#' {
#' 	# Using identify mode
#' 	# Report the data format of the raster file by using the identify mode 
#'  # and specifying a data file name:
#' 	src_dataset <- system.file("external/tahoe_highrez.tif", package="gdalUtils")
#' 	gdalmanage(mode="identify",datasetname=src_dataset)
#' 	
#' 	# Recursive mode will scan subfolders and report the data format:	
#' 	src_dir <- system.file("external/", package="gdalUtils")
#' 	gdalmanage(mode="identify",datasetname=src_dir,r=TRUE)
#' 	
#' \dontrun{
#' 		# Using copy mode	
#' 		# Copy the raster data:
#' 		file_copy <- tempfile(fileext=".tif")
#' 		gdalmanage(mode="copy",src_dataset,file_copy)	
#' 		file.exists(file_copy)
#' 		
#' 		# Using rename mode
#' 		# Rename the raster data:
#' 		new_name <- tempfile(fileext=".tif")
#' 		gdalmanage(mode="rename",file_copy,new_name)	
#' 		file.exists(new_name)
#' 		
#' 		# Using delete mode
#' 		# Delete the raster data:
#' 		gdalmanage(mode="delete",new_name)	
#' 		file.exists(new_name)		
#' }
#' }
#' @export

gdalmanage <- function(
		mode,
		datasetname,newdatasetname,
		r,u,f,
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
					varnames <- c("r","u")),
			vector = list(
					varnames <- c()),
			scalar = list(
					varnames <- c()),
			character = list(
					varnames <- c("mode","f","datasetname","newdatasetname")),
			repeatable = list(
					varnames <- c())
	)
	
	parameter_order <- c(
			"mode","r","u","f","datasetname","newdatasetname")
	
	parameter_noflags <- c("mode","datasetname","newdatasetname")
	
	parameter_doubledash <- NULL
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdalmanage"
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
	
	if(mode=="identify") return(cmd_output) else return(NULL)
}
