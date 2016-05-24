#' Make a matrix symmetric
#' 
#' Make the matrix symmetric by making all "mirrored" positions consistent. A variety of methods are provided to make the matrix symmetrical.
#' @param matrix The matrix to make symmatric
#' @param method The method to use to make the matrix symmetric. Default is to take the maximum. 
#' \itemize{
#' 	\item{"max"} {For each position, \eqn{m_{i,j}}, use the maxiumum of \eqn{(m_{i,j}, m_{j,i})}}
#' 	\item{"min"} {For each position, \eqn{m_{i,j}}, use the minimum of \eqn{(m_{i,j}, m_{j,i})}}
#' 	\item{"avg"} {For each position, \eqn{m_{i,j}}, use the mean: \eqn{(m_{i,j} + m_{j,i})/2}}
#' 	\item{"ld"} {Copy the lower triangular portion of the matrix to the upper triangular portion.}
#' 	\item{"ud"} {Copy the upper triangular portion of the matrix to the lower triangular portion.}
#' }
#' @param adjacencyList Logical. If false, returns the symmetric matrix (the same format as the input). If true, returns an adjacency list representing the upper triangular portion of the adjacency matrix with addressing based on the row.names of the matrix provided.
#' @return The symmetric matrix
#' @description This function was part of the package ENA. This package is no more available
#' and it cannot be installed from archive because some dependencies are no more
#' available.
#' @export
#' @author Jeffrey D. Allen \email{Jeffrey.Allen@@UTSouthwestern.edu}
#' @examples
#' #Create a sample 3x3 matrix
#' mat <- matrix(1:9, ncol=3)
#' 
#' #Copy the upper diagonal portion to the lower
#' symmetricize(mat, "ud")
#' 
#' #Take the average of each symmetric location
#' symmetricize(mat, "avg")
#' 
symmetricize <-
function(matrix, method=c("max", "min","avg", "ld", "ud"), adjacencyList=FALSE){
	method <- match.arg(method)
	
	if (missing(matrix)){
		stop("You must provide a matrix to symmetricize.")
	}
	
	x <- matrix
	if (method=="ld"){
		#solution from Michael Conklin, http://www.biostat.wustl.edu/archives/html/s-news/2000-03/msg00127.html
		x[matrix(c(col(x)[lower.tri(x)], row(x)[lower.tri(x)]), ncol = 2)] <-
			x[matrix(c(row(x)[lower.tri(x)], col(x)[lower.tri(x)]), ncol = 2)]
	}
	if (method=="ud"){				
		x[matrix(c(col(x)[upper.tri(x)], row(x)[upper.tri(x)]), ncol = 2)] <-
			x[matrix(c(row(x)[upper.tri(x)], col(x)[upper.tri(x)]), ncol = 2)]
	}
	if (method=="max" || method=="min" || method=="avg"){
		#retrieve the (row,col) indices for the lower-diagonal entries in the matrix
		ldi <- matrix(c(row(x)[lower.tri(x)], col(x)[lower.tri(x)]), ncol = 2)
			#already column-major, so we can leave this
		
		#retrieve the (row,col) indices for the upper-diagonal entries in the matrix
		udi <- matrix(c(row(x)[upper.tri(x)], col(x)[upper.tri(x)]), ncol = 2)
			#in order for these to be in the symmetrical order as ldi, we need to sort.
		udi <- udi[order(udi[,1], udi[,2]),]
		
		#extract the upper and lower diagonal elements in a way that's symmetrical in their indexces
		ud <- x[udi]
		ld <- x[ldi]
		
		#replace with either the min, max, or mean, depending on the selection
		if (method=="max"){			
			x[ldi] <- apply(rbind(ld,ud),2,max);
			x[udi] <- apply(rbind(ld,ud),2,max);
		}
		if (method=="min"){
			x[ldi] <- apply(rbind(ld,ud),2,min);
			x[udi] <- apply(rbind(ld,ud),2,min);
		}
		if (method=="avg"){
			x[ldi] <- apply(rbind(ld,ud),2,mean);
			x[udi] <- apply(rbind(ld,ud),2,mean);
		}
	}	
	
	if (!adjacencyList){
		#return the adjacency matrix
		return(x)
	}
	else{
		#convert to adjacency list and add addressing
		if (is.null(rownames(matrix))){
			stop("You requested adjacency list format, but the matrix provided has no (row) names.")
		}		
		names <- rownames(matrix)		
		add <- getFromNamespace(".getTableAddressing", ns="HelpersMG")(names)
		toReturn <- (cbind(add, x[upper.tri(x)]))		
		colnames(toReturn) <- c("Source", "Dest", "Symmetric")
		return(toReturn)
	}
}
