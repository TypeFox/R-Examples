NULL
#'
#' ks.test S3 method for 'climdex.data.frame'
#'
#' @param data a \code{\link{climdex.data.frame}} object 
#' @param observed name (String) of the column of \code{data} containing the obseved climate indices
#' @param generated names (String vector) of the columns of \code{data} containing the climate index realizations which will be tested.
#' @param ... further arguments
#' @export

#' @title Kolgomorov-Smirnov Tests for a ClimDex Data Frame
#' @seealso \code{\link{climdex.data.frame}},\code{\link{wilcox.test}},\code{\link{ks.test}}
#' 
#' @author Annalisa Di Piazza, Emanuele Cordano
#'
#' @examples
#' 
#' # See the example of 'climdex.data.frame' function
#'  
#'


ks.test.climdex.data.frame <- function(data,observed,generated,...) {
	
	data <- as.data.frame(data)
	out <- NULL
#	
#	out <- list()
	
	
#	for (it in generated) {
#		
#		index <- paste(observed,"vs",it,sep="_")
#		out[[index]] <- ks.test(x=data[,observed],y=data[,it],...)
#	}
	index <- paste(observed,"vs",generated,sep="_")
	out <- lapply(X=generated,FUN=function(it) {
				ks.test(x=data[,observed],y=data[,it],...)
			
	})

	names(out) <- index
	
	
	return(out)
}
