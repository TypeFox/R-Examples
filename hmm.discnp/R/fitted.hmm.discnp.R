fitted.hmm.discnp <- function(object,...) {
y <- object$y
if(is.null(y)) stop("Observations \"y\" were not kept.\n")
if(!is.list(y)) y <- list(y)
if(!is.numeric(y[[1]]))
if(!all(sapply(y,is.numeric)))
	stop(paste("Some observations are not numeric;\n",
                   "fitted values make no sense.\n"))
sp(y,object,means=TRUE)$means
}
