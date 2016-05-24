#' Get a block of raster cell values (optimization fix for RasterStacks)
#' 
#' A faster version of getValuesBlock for RasterStack.
#' 
#' @param x	Raster* object
#' @param row positive integer. Row number to start from, should be between 1 and nrow(x)
#' @param nrows postive integer. How many rows? Default is 1
#' @param col postive integer. Column number to start from, should be between 1 and ncol(x)
#' @param ncols	postive integer. How many columns? Default is the number of colums left after the start column
#' @param lyrs integer (vector). Which layers? Default is all layers (1:nlayers(x))
#' 
#' @return matrix or vector (if (x=RasterLayer), unless format='matrix')
#' @author Jonathan A. Greenberg
#' @seealso \code{\link[raster]{getValuesBlock}}
#' @details In certain cases, getValuesBlock may run very slowly on a RasterStack,
#' particularly when the RasterStack is comprised of RasterBricks.  This code attempts
#' to fix the inefficiency by running the extract on each unique file of the RasterStack,
#' rather than each unique layer.
#' 
#' @examples
#' tahoe_highrez <- brick(system.file("external/tahoe_highrez.tif", package="spatial.tools"))
#' tahoe_highrez_stack <- stack(tahoe_highrez,tahoe_highrez,tahoe_highrez)
#' # getValuesBlock stack extraction:
#' system.time(tahoe_highrez_extract <- getValuesBlock(tahoe_highrez_stack))
#' # getValuesBlock_stackfix stack extraction:
#' system.time(tahoe_highrez_extract <- getValuesBlock_stackfix(tahoe_highrez_stack))
#' @import raster
#' @export

getValuesBlock_stackfix <- function(x, row=1, nrows=1, col=1, ncols=(ncol(x)-col+1), lyrs=(1:nlayers(x)))
{
	single_filename <- NULL
	if(class(x)=="RasterStack")
	{
		# First we will determine the unique files
		all_filenames <- sapply(x@layers,function(X) { filename(X) } )
		inMemory_layers <- (1:nlayers(x))[sapply(x@layers,function(X) { inMemory(X) } )]
		
		unique_filenames <- unique(all_filenames)
		unique_getValuesBlock <- 
				foreach(single_filename=unique_filenames,.packages=c("raster")) %dopar%
				{
					if(single_filename!="")
						getValuesBlock(brick(single_filename),row,nrows,col,ncols)
					else
						getValuesBlock(stack(x,bands=inMemory_layers),row,nrows,col,ncols)
				}
		
		band_layers <- sapply(x@layers,function(x) x@data@band)
		nlyrs_out <- length(lyrs)
		# This could be sped up a bit with more clever vectorizing but...
		out_matrix <- matrix(nrow=(nrows*ncols),ncol=nlyrs_out)
		for(i in 1:nlyrs_out)
		{
			current_layer <- lyrs[i]
			file_index <- which(all_filenames[i] == unique_filenames)
			out_matrix[,i] <- unique_getValuesBlock[[file_index]][,band_layers[i]]
		}
		return(out_matrix)
	} else
	{
		if(class(x)=="RasterLayer" || nlayers(x)==1)
			return(getValuesBlock(x, row=row, nrows=nrows, col=col, ncols=ncols))
		else
			return(getValuesBlock(x, row=row, nrows=nrows, col=col, ncols=ncols, lyrs=lyrs))
	}
}
