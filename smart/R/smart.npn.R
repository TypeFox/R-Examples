smart.npn <-
function(x, npn.func = "kendall", npn.thresh = NULL, verbose = TRUE){
	gcinfo(FALSE)
	n = nrow(x)
  	d = ncol(x)
  	x.col = colnames(x)
  	x.row = rownames(x)
		y = qnorm(apply(x,2,rank)/(n+1))
		y = y/sd(y[,1])

 
    # identity
  if(npn.func == "pearson"){
  	  x = cor(x)
   		colnames(x) = x.col
		  rownames(x) = x.col
  }  	
  	# Shrinkaage transformation
	if(npn.func == "ns"){
		if(verbose) cat("Conducting the nonparanormal (npn) transformation via normal score....")
		
		x=cor(y)
		
		if(verbose) cat("done.\n")
		rm(n,d,verbose)
   		gc()	
   		colnames(x) = x.col
	  	rownames(x) = x.col
	}
	
	# Truncation transformation
	if(npn.func == "npn"){
		if(verbose) cat("Conducting nonparanormal (npn) transformation via npn....")
		if(is.null(npn.thresh)) npn.thresh = 1/(4*(n^0.25)*sqrt(pi*log(n)))
		
		 x = qnorm(pmin(pmax(apply(x,2,rank)/n, npn.thresh), 1-npn.thresh))
 	   x = x/sd(x[,1])
     x=cor(x)
     	
    	if(verbose) cat("done.\n")
    	rm(n,d,npn.thresh,verbose)
   		gc()
   		colnames(x) = x.col
		  rownames(x) = x.col
	}
	
	if(npn.func == "spearman"){
		if(verbose) cat("Conducting nonparanormal (npn) transformation via spearman....")
		x = 2*sin(pi/6*cor(x,method="spearman"))
    	if(verbose) cat("done.\n")
    	rm(n,d,verbose)
   		gc()
   		colnames(x) = x.col
		  rownames(x) = x.col
	}

	if(npn.func == "kendall"){
		if(verbose) cat("Conducting nonparanormal (npn) transformation via kendall....")
		x = sin(pi/2*cor.fk(x))
    	if(verbose) cat("done.\n")
    	rm(n,d,verbose)
   		gc()
   		colnames(x) = x.col
		  rownames(x) = x.col
	}
	
	return(list(cov=x,scaled=y))
}
