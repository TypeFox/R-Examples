#' Writes image data to a flat binary file using col/row/band positioning.
#' 
#' @param filename Character.  The path and filename of a "blank" binary file to store the image data.
#' @param mode The mode of data on disk.  Defaults to real64() (double precision floating point).
#' @param image_dims Vector. Vector of length(image_dims)==3 representing the number of columns, rows and bands in the output image.
#' @param interleave Character. The require output interleave.  By default is "BSQ". OTHER INTERLEAVES CURRENTLY UNSUPPORTED.
#' @param data Vector, matrix, array, or other data source which is coercible to a vector. This is the data to be written to the image.
#' @param data_position List. A length==3 list, containing the column, row, and band positions ranges to write the output data.
#' 
#' @return NULL
#' @author Jonathan A. Greenberg
#' @seealso \code{\link[mmap]{mmap}},\code{\link[spatial.tools]{create_blank_raster}}
#' @keywords mmap
#' @examples \dontrun{ 
#' tahoe_highrez <- brick(system.file("external/tahoe_highrez.tif", package="spatial.tools"))
#' # Create a blank file using create_blank_raster
#' test_blank_file <- create_blank_raster(reference_raster=tahoe_highrez)
#' blank_raster <- brick(test_blank_file)
#' # It should be all 0s:
#' setMinMax(blank_raster)
#' # Write some ones to to the 100th line, columns 25 to 50, bands 1 and 3:
#' data_position <- list(25:50,100,c(1,3))
#' data1s <- array(1,dim=c(
#' 	length(data_position[[1]]),
#' 	length(data_position[[2]]),
#' 	length(data_position[[3]])))
#' plot(raster(test_blank_file,layer=1))
#' binary_image_write(filename=test_blank_file,
#' 	mode=real64(),image_dims=dim(tahoe_highrez),interleave="BSQ",
#' 	data=data1s,data_position=data_position)
#' setMinMax(blank_raster)
#' plot(raster(blank_raster,layer=1))
#' }
#' @import mmap
#' @export

binary_image_write=function(filename,mode=real64(),image_dims,interleave="BSQ",
	data,data_position)
{
	
	if(is.character(mode))
	{
		mode = 
				#spatial.tools:::
		dataType_converter(from=mode)
	}
	
	# data_position should be of format rbind(col_pos,row_pos,band_pos)
	# UL corner is 1,1,1
	
	if(class(data_position)=="list")
	{
		data_position=
			t(expand.grid(data_position[[1]],data_position[[2]],data_position[[3]]))
	}
	
	if(dim(data_position)[1]==2)
	{
		data_position <- rbind(data_position,rep(1,dim(data_position)[2]))
	}
	
	# Some error checking up here.
	
	image_x=image_dims[1]
	image_y=image_dims[2]
	if(length(image_dims)==3)
	{
		image_z=image_dims[3]
	} else
	{
		image_z=1
	}
	
	if(interleave=="BSQ")
	{
		cell_position=as.integer(
			((data_position[2,]-1)*image_x)+
			(data_position[1,])+
			((data_position[3,]-1)*(image_x*image_y))
			)
	}
	
	if(class(data)=="array")
	{
		data=as.matrix(data,nrow=image_x*image_y,ncol=image_z)
	} 
		
	out = mmap(filename, mode=mode)	
	out[cell_position] <- as.numeric(data)
	munmap(out)	
}