## from the original random forest package
importance <- function(x, ...)  UseMethod("importance")

importance.default <- function(x, ...)
    stop("No method implemented for this class of object")

importance.obliqueRF <- function(x, ...) {
    if (!inherits(x, "obliqueRF"))
        stop("x is not of class obliqueRF")
    if(is.null(x$importance))
    {
    	cat("the importance was not calculated, please rerun the calculation\n");
    	print(x$call);
    	stop("");
    }
    	
    x$importance
}
