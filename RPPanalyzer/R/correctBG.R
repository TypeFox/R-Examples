`correctBG` <-
function(x, method="normexp") {
      elements <- names(x)
      

    if(method=="addmin") {

	## for- and backgorund are vectors
	if(is.vector(x[[1]]) && is.vector(x[[2]])) {

	    ## test if both vectors have the same length
	    if (length(x[[1]]) != length(x[[2]])){
		stop("Expression and background vector must have the same length!")
	    }
	    
	    ## substract background
	    res <- x[[1]] - x[[2]]

	    ## make sure everythinh is positive 
	    res <- ifelse(min(res)<=0, res+abs(min(res))+1, res)

       ## replace expression matrix with bg corrected exp matrix
       x[[1]] <- res
       
	    ## return the result
	    return(x)
	    
	}
	## for- and background are matrices
	else if(is.matrix(x[[1]]) && is.matrix(x[[2]])) {

	    ## dimension of the matrices have to be equal
	    if(any(dim(x[[1]]) != dim(x[[2]]))) {
		stop("Expression and background matrix must have the same dimension!")
	    }

	    ## substract background
	    res <- x[[1]] - x[[2]]

	    ## add min columnwise if min is negative or equal zero
	    res <- apply(res, 2, function(xx) {
		    
		    if(min(xx)<=0)
			return(xx+abs(min(xx))+1)
		    else
			return(xx)		
		})
       ## replace expression matrix with bg corrected exp matrix
       x[[1]] <- res
       
	    ## return result
	    return(x)

	}
	else {
	    stop("Both, expression and backgorund have te be either a vector or a matrix!")
	
	}
    }
    ## limma bg correction method
    else {

	## construct a RG List with both channels have the same values
	RG <- new("RGList", list(R=x[[1]], Rb=x[[2]], G=x[[1]], Gb=x[[2]]))

	## use backgroundCorrect from the limma package
	## passing the method
	RGb <- backgroundCorrect(RG, method=method)

	## use the green channel of the resulting RG list as result, following a suggestion of 
	## Gordon Smyth
	## replace expression matrix with bg corrected exp matrix
       x[[1]] <- RGb$G
	return(x)
    }
}

