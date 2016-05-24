#' Easier-to-use function for grabbing a block of data out of a Raster*.
#' 
#' @param x Raster* Some input Raster* object.
#' @param r1 Numeric. The start row of the chunk.
#' @param r2 Numeric. The end row of the chunk.
#' @param c1 Numeric. The start column of the chunk.
#' @param c2 Numeric. The end row of the chunk.
#' @param lyrs Numeric. Vector of layer IDs.  Defaults to all layers (1:nlayers(x)).
#' @param format Character. See Details. 
#' @param ... Other parameters.
#' 
#' @details This allows for a larger number of output formats to be generated
#' when extracting chunks of data from a Raster* object.  If format="array" (default), 
#' the chunk will be returned in a 3-d array with dimensions representing column,row,and layer.  
#' If "raster", the chunk will be returned as a Raster* object.  If "data.frame", it will
#' be returned as a data.frame.  If "data.frame.dims", it will return a list, where the first
#' component (named "values") is the same as the data.frame when using format="data.frame", and
#' the second component (named "dim") is the dimensions of the extracted chunk.
#' 
#' @return An array or raster object.
#' @author Jonathan A. Greenberg
#' @seealso \code{\link[raster]{getValues}}
#' @examples
#' tahoe_highrez <- brick(system.file("external/tahoe_highrez.tif", package="spatial.tools"))
#' mychunk <- getValuesBlock_enhanced(tahoe_highrez,r1=100,r2=110,c1=20,c2=50)
#' class(mychunk)
#' dim(mychunk)
#' mychunk_raster <- getValuesBlock_enhanced(tahoe_highrez,r1=100,r2=110,c1=20,c2=50,format="raster")
#' mychunk_raster
#' @import raster
#' @export

getValuesBlock_enhanced=function(x,r1=1,r2=nrow(x),c1=1,c2=ncol(x),lyrs=seq(nlayers(x)),format="array",...)
{	
	if(format=="array")
	{
		layer_names=names(x)
		getvalues_raw <- as.numeric(getValuesBlock_stackfix(x,row=r1,nrows=(r2-r1+1),col=c1,ncols=(c2-c1+1),lyrs=lyrs))		
		getvalues_raw_nrows=r2-r1+1
		getvalues_raw_ncols=c2-c1+1
		getvalues_raw_nlayers=nlayers(x)
		
		# Test the input file.
		if(getvalues_raw_nlayers==1)
		{
			# Raster
			getvalues_array=array(data=getvalues_raw,
					dim=c(getvalues_raw_ncols,getvalues_raw_nrows,getvalues_raw_nlayers))		
		} else
		{
			# Brick or stack
			getvalues_array=array(data=getvalues_raw,
					dim=c(getvalues_raw_ncols,getvalues_raw_nrows,getvalues_raw_nlayers))
		}
		dimnames(getvalues_array) <- list(NULL,NULL,NULL)
		if(!is.null(layer_names)) dimnames(getvalues_array)[[3]]=layer_names
		return(getvalues_array)
	}
	
	if(format=="raster")
	{
		return(crop(x, extent(x, r1=r1, r2=r2, c1=c1,c2=c2)))
	}
	
	if(format=="data.frame" || format=="data.frame.dims")
	{
		#	layer_names=names(x)
		getvalues_raw <- getValuesBlock_stackfix(x,row=r1,nrows=(r2-r1+1),col=c1,ncols=(c2-c1+1),lyrs=lyrs)	
		
		getvalues_df <- as.data.frame(getvalues_raw)
		names(getvalues_df) <- names(x)
		
		# Fix factors:
		factor_layers <- is.factor(x)
		if(any(factor_layers))
		{
			factor_levels <- levels(x)
			for(i in seq(nlayers(x))[factor_layers])
			{
				temp_factor_levels <- factor_levels[[i]][[1]]
				temp_factor_levels <- data.frame(ID=temp_factor_levels[,1],code=temp_factor_levels[,2])
				temp_df_data <- getvalues_df[[i]]
				
				temp_factor_column <- with(temp_factor_levels,code[match(temp_df_data,ID)])
				
				getvalues_df[[i]] <- temp_factor_column
			}
		}
		
		names(getvalues_df) <- names(x)
		
		if(format=="data.frame.dims")
		{
			getvalues_raw_nrows=r2-r1+1
			getvalues_raw_ncols=c2-c1+1
			getvalues_raw_nlayers=nlayers(x)
			getValues.dims <- c(getvalues_raw_ncols,getvalues_raw_nrows,getvalues_raw_nlayers)
			
			getvalues_df <- list(values=getvalues_df,dim=getValues.dims)
		}
	
		return(getvalues_df)
		
	}
}