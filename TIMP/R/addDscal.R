"addDscal" <-
function (modellist, dscal) 
{
    #dscal is list with structure:
    # list ( list (to, from, value = c(vector of scaling),perclp =logical))

	if(length(dscal) != 0) {
	  for(i in 1:length(dscal)){
	      if(length(dscal[[i]]$perclp) == 0)
		dscal[[i]]$perclp <- FALSE
	      modellist[[dscal[[i]]$to]]@dscalspec <- dscal[[i]] 
	      modellist[[dscal[[i]]$to]]@drel <- dscal[[i]]$value
	  }
	}
	modellist
}

