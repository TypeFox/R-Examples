#' 'Not Available' / Missing Values for ff vectors
#'
#' The generic function \code{is.na} indicates which elements are missing.\cr
#' The generic function \code{is.na<-} sets elements to \code{NA}. 
#'
#' @rdname is.na.ff
#' @method is.na ff 
#' @method is.na<- ff
#' @usage \method{is.na}{ff} (x, ...)
#' @example ../examples/isNA.R
#' @param x a \code{ff} vector
#' @param ... other parameters passed on to chunk
#' @param value a suitable ff index vector for use with x
#' @return A logical \code{ff} vector of the same length of x indicating if the ff vector contains missing values. 
#' @export
#' @export is.na.ff
#' @seealso \code{\link[base]{is.na}, \link[ff]{ffvecapply}}
is.na.ff <- function(x, ...){
	res <- ff(vmode="logical", length=length(x))
	for (i in chunk(x, ...)){		
		res[i] <- is.na(x[i])		
	}
	res		
}

#' @rdname is.na.ff
#' @usage \method{is.na}{ff} (x, ...) <- value
#' @export
"is.na<-.ff" <- function(x, ..., value){
	if(inherits(value, "ff_vector")){
		for (i in chunk(value, ...)){
  		set.ff(x=x, i=value[i], value=NA, add = FALSE)
  	}	
	}else{
		set.ff(x=x, i=value, value=NA, add = FALSE)
	}	
	x	
}


