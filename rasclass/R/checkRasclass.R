##################################################################################
# Set Method: checkRasclass
##################################################################################
checkRasclass <- function(object){}
setMethod('checkRasclass', signature(object = 'rasclass'),

function(object){
	
	allOk <- TRUE
	
	# Check for skeleton consistency
	if(is.null(object@gridSkeleton@nrows)){
		cat('Warning: Number of rows in gridSkeleton is not specified\n')
		allOk <- F
	}
	else if(is.null(object@gridSkeleton@ncols)){
		cat('Warning: Number of columns in gridSkeleton is not specified\n')
		allOk <- FALSE
	}
	else if(object@gridSkeleton@nrows*object@gridSkeleton@ncols != length(object@gridSkeleton@grid)){
		cat('Warning: The number of cells does not correspond with the number of rows and cols in gridSkeleton.\n')
		allOk <- FALSE
	}

	# Check if skeleton is complete
	if(is.null(object@gridSkeleton@xllcorner)){
		cat('Warning: xllcorner in gridSkeleton is not specified\n')
		allOk <- FALSE	
	}
	if(is.null(object@gridSkeleton@yllcorner)){
		cat('Warning: yllcorner in gridSkeleton is not specified\n')
		allOk <- FALSE	
	}
	if(is.null(object@gridSkeleton@cellsize)){
		cat('Warning: cellsize in gridSkeleton is not specified\n')
		allOk <- FALSE	
	}
	if(is.null(object@gridSkeleton@NAvalue)){
		cat('Warning: NAvalue in gridSkeleton is not specified\n')
		allOk <- FALSE	
	}

	# Ceck for consistency of skeleton vs data
	if(is.null(object@gridSkeleton@grid)){
		cat('Warning: Indicator grid in gridSkeleton is not specified\n')
		allOk <- FALSE			
	}
	else if(sum(object@gridSkeleton@grid) != nrow(object@data)){
		cat('Warning: The number of observations in prediction data is not consistent with the gridSkeleton\n')
		allOk <- FALSE
	}

	# Check if samplename exists in Data
	if(!is.element(object@samplename, names(object@data))){
		cat('Warning: Samplename specified wrongly, "', object@samplename, '" is not a column in the rassclass object data\n', sep = '')
		allOk <- FALSE
	}
	
	# Return value
	allOk
}
)