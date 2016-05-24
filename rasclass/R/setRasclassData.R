##################################################################################
# Set Method: setRasclassData
##################################################################################
setRasclassData <- function(newdata, object = new('rasclass'), ncols = NA, nrows = NA, xllcorner = NA, yllcorner = NA, cellsize = NA, NAvalue = NA, samplename = 'sample'){}

setMethod('setRasclassData', signature(newdata = 'data.frame'),

function(newdata, object = new('rasclass'), ncols = NA, nrows = NA, xllcorner = NA, yllcorner = NA, cellsize = NA, NAvalue = NA, samplename = 'sample'){
	
	# Set path
	object@path <- 'Data specified manually using setRasclassData()'
	
	# Set sample name
	object@samplename <- samplename
	
	# Remove data path
	object@path <- as.character(NA)
	
	# Update the gridSkeleton
	if(!is.na(NAvalue)) 	object@gridSkeleton@NAvalue   <- NAvalue
	if(!is.na(ncols)) 		object@gridSkeleton@ncols     <- ncols
	if(!is.na(nrows)) 		object@gridSkeleton@nrows     <- nrows
	if(!is.na(xllcorner)) 	object@gridSkeleton@xllcorner <- xllcorner
	if(!is.na(yllcorner)) 	object@gridSkeleton@yllcorner <- yllcorner
	if(!is.na(cellsize)) 	object@gridSkeleton@cellsize  <- cellsize
	
	# Set up the gridSkeleton grid (NAhandle)
	object@gridSkeleton@grid <- as.integer(unlist(apply(newdata[names(newdata) != samplename], 1, function(x) !is.element(NA, x))))
	
	# Remove NAs and set data
	object@data <- newdata[as.logical(object@gridSkeleton@grid), ]
	
	# Convert sample to factor
	object@data[, object@samplename] <- factor(object@data[, object@samplename])
	
	# Build Formula
	object <- buildFormula(object)
	
	# Check consistency
	if(!checkRasclass(object)){
		stop('Data object is not consistent, check data and try again')
	}
	
	# Return object
	object
}
)