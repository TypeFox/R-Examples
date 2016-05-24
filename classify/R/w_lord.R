
wlord <- function(probs=NULL, cats=NULL){

	#Wrapper for Wingersky Lord
	if(!is.numeric(cats)){
		stop("Item categories need to be specified as a numeric vector.")		
	}
	if(!is.matrix(probs)){
		stop("Item probabilities need to be specified as a matrix.")
	}
	if(nrow(probs)!=length(cats)){
		stop(paste("Categories Specified for:",length(cats), "items, but probabilities provided for" ,nrow(probs),"items."))
	}
	if(min(cats)<2){
		stop("Minimum number of item categories is 2.")
	}	
	
	ret <- .Call( "rcpp_w_lord_c", probs, cats,PACKAGE = "classify" )
	return(ret)
}