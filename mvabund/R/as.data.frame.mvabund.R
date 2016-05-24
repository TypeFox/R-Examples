###############################################################################
# as.data.frame.mvabund: converge a mvabund object to a data.frame     		    #
# this function should be invisible, can only be used through as.data.frame   #
###############################################################################

as.data.frame.mvabund <- function(x, row.names = NULL, optional = FALSE, ...)
{
	if(!is.mvabund(x))
    stop("The function 'as.data.frame.mvabund' can only be used for a mvabund object.")

  if(is.null(dim(x))) x <- mvabund(c(x))

  classx <- class(x)
	
	mc <- match.call(expand.dots=FALSE)
	
	dots <- mc$...
	
	if (length(dots)>0){
     mc$... <- NULL
     dots <- lapply( dots, eval, parent.frame() )
  }
  
  if(length(dim(x))>2) {
        class(x) <- "array"
        return(as.data.frame(x))
        mceval <- lapply(mc , eval, parent.frame())
        # Delete the function body in the call.
	      mceval[[1]] <- NULL
	      mceval$x    <- NULL
	      # Converge x into a data.frame.
	      x <- do.call( "as.data.frame", c(list(x), mceval, dots))
        return(x)
  }
	
	if(length(classx)==1) {
	# x is has no other class than mvabund.
		n <- NROW(x)
		p <- NCOL(x)
		# Converge x into a matrix so that it can be further converged
    # to a data frame.
		x <- matrix(x, nrow=n, ncol=p)
		} 
	# Use either the created matrix class or the second class of x to
  # find a coercion method.
	x <- unabund(x)

	mceval <- lapply(mc , eval, parent.frame())
	
  # Delete the function body in the call.
	mceval[[1]] <- NULL
	mceval$x    <- NULL
	# Converge x into a data.frame.
	x <- do.call( "as.data.frame", c(list(x), mceval, dots))
	
	return(x)
	
}

