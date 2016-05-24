`cem.summary` <-
function (obj, verbose = 0) 
{
    if (is.null(obj$matched)) 
	return(NULL)
    tab <- matrix(, 3, obj$n.groups)
    tab[1, ] <- table(obj$groups)
    id1 <- numeric(obj$n.groups)
    id0 <- numeric(obj$n.groups)
    for(i in 1:obj$n.groups){
		id1[i] <- length(which(obj$groups==obj$g.names[i] & obj$matched==TRUE))
		id0[i] <- length(which(obj$groups==obj$g.names[i] & obj$matched==FALSE))
    }
	
    tab[2,] <- id1
    tab[3,] <- id0
	
	
    colnames(tab) <- paste("G", obj$g.names, sep = "")
    rownames(tab) <- c("All", "Matched", "Unmatched")
    if (verbose > 1) {
        cat(sprintf("\nCEM Subclasses: %d\n", length(obj$mstrataID)))
        cat("\nSample sizes:\n")
        print(tab)
    }
    return(tab)
}