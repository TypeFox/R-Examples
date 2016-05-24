msgpack.matrix <-
function(data) {
	n <- 1
	n_row <- length(data)
	n_col <- length(data[[1]])
	nms <- names(data[[1]])
	result <- list()
	for ( row in data ) {
		for ( col in row ) {
			result[n] <- col
			n <- n+1
		}
	}

	if ( n_col == 1 && n_row == 1 ) {
		result <- result[[1]]
		if ( !is.null(nms) ) {
			names(result) <- nms
		}
	}
	else if ( n_row > 1 ) {
		result <- matrix(result, n_row, n_col, byrow=TRUE)
		if ( !is.null(nms) ) {
			colnames(result) <- nms
		}
	}
	
	return(result)
}
