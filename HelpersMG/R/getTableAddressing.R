#' Get the adjacency list addressing template. 
#' 
#' Useful if you want to store the networks in their condensed upper-diagonal form while still
#' having the benefit of convenient addressing and/or if you are using a simulated dataset in
#' which you know the truth and want to store all the values in a single data.frame.
#' 
#' Internal function used to get the addressing template for a data.frame to contain the adjacency 
#' list representation of a matrix.
#' @param variableNames the names of all genes to include in the adjacency list
#' @param truth The true adjacency matrix. Often will not be available, but is useful for 
#' debugging and testing.
#' @return A data.frame representing the adjacency list of the matrix provided.
#' @description This function was part of the package ENA. This package is no more available
#' and it cannot be installed from archive because some dependencies are no more
#' available.
#' @export
#' @author Jeffrey D. Allen \email{Jeffrey.Allen@@UTSouthwestern.edu}
#' @examples
#' # getTableAddressing example
#' 
.getTableAddressing <- function (variableNames, truth){		
	joint <- matrix(0, nrow=length(variableNames), ncol=length(variableNames))
	rownames(joint) <- variableNames
	colnames(joint) <- variableNames
	
	#extract the row, col indices for the upper-triangular portion of the matrix
	address <- cbind(row(joint)[upper.tri(joint)], col(joint)[upper.tri(joint)])
	address <- matrix(row.names(joint)[address], ncol=2)
	
	upper <- upper.tri(joint);
	
	if (missing(truth)){
		aggregate <- data.frame(address)
		colnames(aggregate) <- c("Source", "Dest")
	}
	else{
		aggregate <- data.frame(address, truth[upper])
		colnames(aggregate) <- c("Source", "Dest", "Truth")
	}
	return(aggregate);
}
