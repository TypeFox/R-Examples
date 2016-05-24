gpcm.rc <- function(beta=NULL,theta=NULL,cats=NULL,alpha=NULL){

	# Do some basic checks
	# Number of categories should equal number of items
	if (is.null(nrow(beta))) {
		#Only one item
		stop("Beta must be a matrix")
	} else {
		if (nrow(beta)!=length(cats)){
			stop("Number of categories should equal number of items")
		}
	}	
	#If alpha is null use PCM
	if (is.null(alpha)){
		alpha <- rep(1,length(cats))
	}
	#Alpha should be same length as cats
	if (length(alpha)!=length(cats)){
		stop("Number of alpha values should equal number of items")
	}
	
	ret <- .Call( "gpcm_c", beta,theta,cats,alpha,PACKAGE = "classify" )
	return(ret)
}

expected.rc <- function(beta=NULL,theta=NULL,cats=NULL,alpha=NULL){

	# Do some basic checks
	# Number of categories should equal number of items
	if (is.null(nrow(beta))) {
		#Only one item
		stop("Beta must be a matrix")
	} else {
		if (nrow(beta)!=length(cats)){
			stop("Number of categories should equal number of items")
		}
	}	
	#If alpha is null use PCM
	if (is.null(alpha)){
		alpha <- rep(1,length(cats))
	}
	#Alpha should be same length as cats
	if (length(alpha)!=length(cats)){
		stop("Number of alpha values should equal number of items")
	}
	
	ret <- .Call( "exp_c", beta,theta,cats,alpha,PACKAGE = "classify" )
	return(ret)
}