betapart.core <- function(x){
	
	# test for a binary numeric matrix or data.frame
	if(! is.matrix(x)){
		x<-as.matrix(x)
  	}
	
	if(! is.numeric(x))
    	stop("The data in x is not numeric.",call.=TRUE)
	
	# simple test for binary data
	xvals <-  unique(as.vector(x))
	if (any(!is.element(xvals, c(0, 1)))) 
        stop("The table contains values other than 0 and 1: data should be presence/absence.", call. = TRUE)

	shared <- x %*% t(x)
    not.shared <-  abs(sweep(shared, 2, diag(shared)))
		
	sumSi <- sum(diag(shared)) # species by site richness
    St <- sum(colSums(x) > 0)  # regional species richness
    a <- sumSi - St            # multi site shared species term

    sum.not.shared <- not.shared + t(not.shared)
    max.not.shared <- pmax(not.shared, t(not.shared))
    min.not.shared <- pmin(not.shared, t(not.shared))

	computations<-list(data=x, sumSi=sumSi, St=St, a=a, shared=shared, not.shared=not.shared, 
	                   sum.not.shared=sum.not.shared, max.not.shared=max.not.shared, 
	                   min.not.shared=min.not.shared)
    class(computations)<-"betapart"

	return(computations)
} 