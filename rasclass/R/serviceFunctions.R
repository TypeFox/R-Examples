##################################################################################
# Set Method: summary
##################################################################################
setMethod('summary', signature(object = 'rasclass'),

function(object){
	
	cat('Call:\n')
	show(object@call)

	cat('\nFormula:\n')
	cat(object@formula)
	
	cat('\n\nData:\n')
	show(summary(object@data))

	cat('\nOverall Accuracy:', object@overallAccuracy)
	cat('\n\nKappa Coefficient:', object@kappa)
	
	cat('\n\nAccuarcy Table:\n')
	show(object@accuracyMatrix)
}
)

##################################################################################
# Set Method: image
##################################################################################
setMethod('image', signature(x = 'rasclassRaster'),

function(x, ...){	
	mymatrix <- matrix(x@grid, nrow = x@ncols)
	mymatrix <- t(apply(mymatrix, 1, rev))
	X <-  x@xllcorner + x@cellsize * 1:x@ncols # Longitude
	Y <-  x@yllcorner + x@cellsize * 1:x@nrows # Latitude
	image(X, Y, mymatrix, col = rainbow(length(levels(factor(x@grid)))), ...)
}
)

##################################################################################
# Set Method: image
##################################################################################
setMethod('image', signature(x = 'rasclass'), function(x){image(x@predictedGrid)})

##################################################################################
# Set Method: View
##################################################################################
setMethod('View', signature(x = 'rasclass'), function(x, title){View(x@data, title)})
