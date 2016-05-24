#' Returns significance level.
#' 
#' Returns the significance level as stars, or NA if a non-numeric value is
#' passed in.
#' 
#' @param x p-value.
star <- function(x) {
	if(is.numeric(x)) {
		return(cut(x, 
				   breaks=c(-Inf,.001,.01,.05,.1,Inf), 
				   labels=c('***','**','*','.',''))
		)
	} else {
		return(NA)
	}
}
