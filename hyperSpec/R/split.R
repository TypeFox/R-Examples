.split <- function (x, f, drop = TRUE){
	validObject (x)
	
	hyperlist <- split (seq (x, index = TRUE), f, drop)
	
	for (i in seq_along (hyperlist)){
		hyperlist [[i]] <- x [hyperlist [[i]], ]
	}
	
	hyperlist
}

##' Split a hyperSpec object according to groups
##' \code{split} divides the \code{hyperSpec} object into a list of
##' \code{hyperSpec} objects according to the groups given by \code{f}.
##' 
##' The \code{hyperSpec} objects in the list may be bound together again by
##' \code{\link{bind} ("r", list_of_hyperSpec_objects)}.
##' 
##' @name split
##' @rdname split
##' @aliases split split-methods split,ANY-method split,hyperSpec-method
##' @docType methods
##' @param x the \code{hyperSpec} object
##' @param f a factor giving the grouping (or a variable that can be converted
##'   into a factor by \code{as.factor})
##' @param drop if \code{TRUE}, levels of\code{f} that do not occur are
##'   dropped.
##' @return A list of \code{hyperSpec} objects.
##' @author C. Beleites
##' @seealso \code{\link[base]{split}}
##' @keywords methods
##' @export
##' @examples
##' 
##' dist <- pearson.dist (chondro[[]])
##' dend <- hclust (dist, method = "ward")
##' z <- cutree (dend, h = 0.15)
##' 
##' clusters <- split (chondro, z)
##' length (clusters)
##' 
##' # difference in cluster mean spectra
##' plot (apply (clusters[[2]], 2, mean) - apply (clusters[[1]], 2, mean))
##' 
##' 


setMethod ("split", signature = signature (x = "hyperSpec"), .split)
