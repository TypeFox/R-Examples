#' get_subdatasets
#' 
#' Returns HDF4, HDF5, and NetCDF subdataset names for standardized files.
#' 
#' @param datasetname Character. Input HDF4/5 or NetCDF file.
#' @param names_only Logical. Return subdataset names only?  Default=TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  

#' @return character vector of subdataset names that can be used in gdal_translate.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) and Matteo Mattiuzzi (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details Currently, this only returns the subdataset names of HDF4, HDF5, and NetCDF files,
#' assuming they follow the SUBDATASET_n_NAME convention.  This can be used with gdal_translate
#' to extract a single subdataset (or with gdal_translate(...,sd_index=n)
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#'
#' @references \url{http://www.gdal.org/gdalinfo.html}
#' 
#' @examples 
#' \dontrun{ 
#' hdf4_dataset <- system.file("external/test_modis.hdf", package="gdalUtils")
#' get_subdatasets(hdf4_dataset)
#' }
#' @importFrom utils glob2rx
#' @export

get_subdatasets <- function(datasetname,names_only=TRUE,verbose=FALSE)
{
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation()
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	if(names_only)
	{
		gdalinfo_raw <- gdalinfo(datasetname)
		subdataset_rawnames <- gdalinfo_raw[grep(glob2rx("*SUBDATASET*NAME*"),gdalinfo_raw)]
		subdataset_names <- sapply(X=seq(length(subdataset_rawnames)),FUN=function(X)
				{
					split1 <- strsplit(subdataset_rawnames[X],"=")
					return(gsub("\"","",split1[[1]][2]))	
				})
		return(subdataset_names)
	}	
}