#' batch_gdal_translate
#' 
#' Runs gdal_translate on a batch of files
#' 
#' @param infiles Character. A directory or a character vector of files (including their path).  If a directory, all files matching the pattern will be converted.
#' @param outdir Character. Output directory to save the output files.
#' @param outsuffix Character. The suffix to append to the input filename (minus its extension) to generate the output filename(s).
#' @param pattern Character. If infiles is a directory, this is used to limit the file it is searching for.
#' @param recursive Logical. If infiles is a directory, should files be searched for recursively?
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' @param ... Parameters to pass to \code{\link{gdal_translate}}
#' 
#' @return Either a list of NULLs or a list of RasterBricks depending on whether output_Raster is set to TRUE.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net})
#' 
#' @details This function is designed to run gdal_translate in batch mode.  Files are
#' passed to the function either directly as a character vector of filenames, or by 
#' passing it a directory and (typically) a search pattern (e.g. pattern=".tif").  gdal_translate
#' will execute based on parameters passed to it, and the output file will be named based on the
#' input file (stripped of its extension), with the outsuffix appended to it.
#' 
#' If a parallel engine is started and registered with foreach, this program will run in parallel
#' (one gdal_translate per worker).  
#'
#' @references \url{http://www.gdal.org/gdal_translate.html}
#' @seealso \code{\link{gdal_translate}}, \code{\link{list.files}}
#' 
#' @examples \dontrun{ 
#' input_folder <- system.file("external",package="gdalUtils")
#' list.files(input_folder,pattern=".tif")
#' output_folder <- tempdir()
#' # library(spatial.tools)
#' # sfQuickInit() # from package spatial.tools to launch a parallel PSOCK cluster
#' batch_gdal_translate(infiles=input_folder,outdir=output_folder,
#' 	outsuffix="_converted.envi",of="ENVI",pattern=".tif$")
#' list.files(output_folder,pattern="_converted.envi$")
#' # sfQuickStop() # from package spatial.tools to stop a parallel PSOCK cluster
#' }
#' @export
#' @import foreach

# TODO: For nested folders, copy the folder hierarchy?  
# TODO: What to do with duplicate filenames (if recursive=TRUE or if the file stripped of an extension is identical)? 

batch_gdal_translate <- function(infiles,outdir,outsuffix="_conv.tif",pattern=NULL,recursive=FALSE,
		verbose=FALSE,
		...)

{
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation()
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# These are just to avoid errors in CRAN checks.
	infile <- NULL
	outfile <- NULL
	
	if(missing(outdir) || !file.exists(outdir) || !file.info(outdir)$isdir)
	{
		stop("Please select a valid outdir.")	
	}
	
	if(file.info(infiles)$isdir)
	{
		if(verbose) message("Input is a directory...")
		infiles <- list.files(infiles,pattern=pattern,full.names=TRUE,recursive=recursive)
	}
	
	# Check for file existence here...
	checkFiles <- sapply(infiles,function(x) return(file.exists(x)))
	
	if(!all(checkFiles))
	{
		message("The following files are missing:")
		message(infiles[!checkFiles])
		stop()
	}
	
	# Generate output names.
	outfiles <- sapply(infiles,function(x,outsuffix,outdir=outdir)
			{
				outfilename <- paste(remove_file_extension(basename(x)),outsuffix,sep="")
				outfilepath <- normalizePath(file.path(outdir,outfilename),mustWork=FALSE)
			},outsuffix=outsuffix,outdir=outdir)
	
	outputs <- foreach(infile=infiles,outfile=outfiles,.packages="gdalUtils") %dopar%
			{
				gdal_translate(src_dataset=infile,dst_dataset=outfile,verbose=verbose,...)
			}
	
	return(outputs)
}