##################################################################################
# Set Method: readRaster
##################################################################################
readRaster <- function(path, asInteger = FALSE){}
setMethod('readRaster', signature(path = 'character'),

function(path, asInteger = FALSE){
	
	object <- new('rasclassRaster')
	
	temp <- readLines(path, n = 6)

	for(i in 1:6){
		tempstring <- unlist(strsplit(temp[i], split = ' '))
		temp[i] <- tempstring[length(tempstring)]
	}

	object@ncols     <- as.numeric(temp[1])
	object@nrows     <- as.numeric(temp[2])
	object@xllcorner <- as.numeric(temp[3])
	object@yllcorner <- as.numeric(temp[4])
	object@cellsize  <- as.numeric(temp[5])
	object@NAvalue   <- as.numeric(temp[6])
	
	if(asInteger){
		object@grid <- as.integer(round(scan(file = path, skip = 6, na.strings = as.character(object@NAvalue))))
	} else {
		object@grid <- scan(file = path, skip = 6, na.strings = as.character(object@NAvalue))
	}
	
	object
}
)