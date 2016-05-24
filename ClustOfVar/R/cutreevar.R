cutreevar <-
function(obj,k=NULL,matsim=FALSE) {
	#obj as produced by ‘hclustvar’
	#k: an integer scalar 
	cl <- match.call()
	if (is.null(n1 <- nrow(obj$merge)) || n1 < 1) 
		stop("invalid 'tree' (merge component)")
    	n <- n1 + 1
	if (is.null(k)) 
        	stop("'k' must be specified")
	k <- as.integer(k)
	if (k < 1 || k > n) 
            stop(gettextf("elements of 'k' must be between 1 and %d", 
                n), domain = NA)
	part<-cutree(obj,k)
	res <- descript(part,obj$rec,cl,matsim)
	class(res) <- "clustvar"
	return(res)
}

