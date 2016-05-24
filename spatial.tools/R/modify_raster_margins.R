#' Add/subtract rows and columns from Raster*
#' 
#' 
#' @param x A Raster* object.
#' @param extent_delta Numeric vector. How many rows/columns to add/subtract to the left,right,top, and bottom of an image.  Default is c(0,0,0,0) (no change).
#' @param value Value to fill in when adding rows/columns.
#' @return A Raster* object.
#' @author Jonathan A. Greenberg
#' @details A quick way to add/subtract margins from a Raster* object.  extent_delta is a four-element integer vector that 
#' describes how many rows/columns to add to the (left,right,top,bottom) of the image (in that order).  Negative values remove rows,
#' positive values add rows.
#' 
#' @examples
#' tahoe_highrez <- brick(system.file("external/tahoe_highrez.tif", package="spatial.tools"))
#' dim(tahoe_highrez)
#' # Remove one row and column from the top, bottom, left, and right:
#' tahoe_highrez_cropped <- modify_raster_margins(x=tahoe_highrez,extent_delta=c(-1,-1,-1,-1))
#' dim(tahoe_highrez_cropped)
#' # Add two rows to the top and left of the raster, and fill with the value 100.
#' tahoe_highrez_expand <- modify_raster_margins(x=tahoe_highrez,extent_delta=c(2,0,2,0),value=100)
#' dim(tahoe_highrez_expand)
#' @import raster
#' @export

modify_raster_margins <- function(x,extent_delta=c(0,0,0,0),value=NA)
{
	x_extents <- extent(x)
	res_x <- res(x)
	
	x_modified <- x
	
	if(any(extent_delta < 0))
	{
		# Need to crop
		# ul:
		ul_mod <- extent_delta[c(1,3)] * res_x
		ul_mod[ul_mod > 0] <- 0
		lr_mod <- extent_delta[c(2,4)] * res_x
		lr_mod[lr_mod > 0] <- 0
	# This works fine, but for some reason CRAN doesn't like it:	
	#	crop_extent <- as.vector(x_extents)
		crop_extent <- c(x_extents@xmin,x_extents@xmax,x_extents@ymin,x_extents@ymax)
		crop_extent[c(1,3)] <- crop_extent[c(1,3)] - ul_mod
		crop_extent[c(2,4)] <- crop_extent[c(2,4)] + lr_mod
		
		x_modified <- crop(x_modified,crop_extent)
	}
	
	if(any(extent_delta > 0))
	{
		# Need to crop
		# ul:
		ul_mod <- extent_delta[c(1,3)] * res_x
		ul_mod[ul_mod < 0] <- 0
		lr_mod <- extent_delta[c(2,4)] * res_x
		lr_mod[lr_mod < 0] <- 0
#		Again, a hack for CRAN?		
#		extend_extent <- as.vector(x_extents)
		extend_extent <- c(x_extents@xmin,x_extents@xmax,x_extents@ymin,x_extents@ymax)
		extend_extent[c(1,3)] <- extend_extent[c(1,3)] - ul_mod
		extend_extent[c(2,4)] <- extend_extent[c(2,4)] + lr_mod
		
		x_modified <- extend(x_modified,extend_extent,value=value)
	}
	
	return(x_modified)
}