setMethod("nobs", signature(object="mix"),
	function(object, ...) {
		nt <- sum(object@ntimes)
    nmiss <- rep(1,nt)
    #for(i in 1:length(object@response)) {
      for(j in 1:length(object@response[[1]])) {
        nmiss <- nmiss*as.numeric(apply(as.matrix(object@response[[1]][[j]]@y),1,function(x) !any(is.na(x))))
      }
    #}
    n <- sum(nmiss)
		#n <- sum(!apply(object@response[[1]]y,1,function(x) any(is.na(x))))
		if(n!=nt) warning("missing values detected; nobs is number of cases with complete data")
		return(n)
	}
)