##################################################################################
# Set Method: writeRaster
##################################################################################
writeRaster <- function(object, path = 'predictedGrid.asc'){}
setMethod('writeRaster', signature(object = 'rasclassRaster'),

function(object, path = 'predictedGrid.asc'){
	
	temp <- paste('ncols         ', object@ncols)
	temp <- append(temp, paste('nrows         ', object@nrows))
	temp <- append(temp, paste('xllcorner     ', object@xllcorner))
	temp <- append(temp, paste('yllcorner     ', object@yllcorner))
	temp <- append(temp, paste('cellsize      ', object@cellsize))
	temp <- append(temp, paste('NODATA_value  ', object@NAvalue))

	allrows <- split(object@grid, rep(1:object@ncols,each=object@nrows))
	allrows <- sapply(allrows, function(x) {
		x[is.na(x)] <- object@NAvalue
		x <- paste(x, collapse = ' ')
		x
		})

	temp <- c(temp, allrows)
	
	writeLines(temp, con = path)
}
)