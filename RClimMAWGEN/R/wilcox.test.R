NULL
#'
#' wilcox.test S3 method for 'climdex.data.frame'
#' 
#' @param x a \code{\link{climdex.data.frame}} object 
#' @param observed name (String) of the column of \code{data} containing the obseved climate indices
#' @param generated names (String vector) of the columns of \code{data} containing the climate index realizations which will be tested.
#' @param ... further arguments
#' 
#' 
#' @rdname wilcox.test
#' @method wilcox.test climdex.data.frame
#' @S3method wilcox.test climdex.data.frame
#' @aliases wilcox.test wilcox.test.climdex.data.frame 
#' @name wilcox.test
#' @export
#' @import stats
#' @author  Emanuele Cordano, Annalisa Di Piazza
#' 
#' @title Wilcoxon Rank Sum and Signed Rank Tests a ClimDex Data Frame
#' @seealso \code{\link{climdex.data.frame}},\code{\link{ks.test}},\code{\link{ks.test.climdex.data.frame}}
#' @examples
#' 
#' # See the example of 'climdex.data.frame' function
#'  




wilcox.test.climdex.data.frame <- function(x,observed,generated,...) {
	
	out <- NULL
	x <- as.data.frame(x)
#	out <- list()
	
	
#	for (it in generated) {
#		
#		index <- paste(observed,"vs",it,sep="_")
#		out[[index]] <- wilcox.test(x=data[,observed],y=data[,it],...)
#	}
#	

	index <- paste(observed,"vs",generated,sep="_")
	

	out <- lapply(X=generated,FUN=function(it,x,observed,...) {
				
				return(wilcox.test(x=x[,observed],y=x[,it],...))
	
			},x=x,observed=observed,...)
	
	names(out) <- index
	
	return(out)
}
