matrix.x<-function (comm, traits, scale = TRUE, notification = TRUE) 
{
	comm<-as.matrix(comm)
	traits<-as.matrix(traits)
    matrix.w <- sweep(comm, 1, rowSums(comm, na.rm=TRUE), "/")
	w.NA <- apply(matrix.w, 2, is.na)
    matrix.w[w.NA] <-0
    if(notification==TRUE){
    	if(length(which(unique(as.vector(w.NA))==TRUE))>0)
    	{
			warning("Warning: NA in community data",call.=FALSE)		
    	}
    }
    x.NA <- apply(traits, 2, is.na)
    if(notification==TRUE){
    	if(length(which(unique(as.vector(x.NA))==TRUE))>0)
    	{
			warning("Warning: NA in traits matrix",call.=FALSE)		
    	}
    }
    
    if (scale == "TRUE") {
        dist.traits <- vegdist(traits, method = "gower", diag = TRUE, upper = TRUE,na.rm=TRUE)
        similar.traits <- 1 - as.matrix(dist.traits)
        matrix.traits <- 1/colSums(similar.traits,na.rm=TRUE)
        matrix.u <- sweep(similar.traits, 1, matrix.traits, "*")
    }
    else{
    	dist.traits <- as.matrix(vegdist(traits, method = "euclidean", diag = TRUE, upper = TRUE,na.rm =TRUE))
	    similar.traits <- 1 - (dist.traits/max(dist.traits,na.rm=TRUE))
    	matrix.traits <- 1/colSums(similar.traits,na.rm=TRUE)
	    matrix.u <- sweep(similar.traits, 1, matrix.traits, "*")
    }
	u.NA <- apply(matrix.u, 2, is.na)
    if (notification == TRUE) {
        if (length(which(unique(as.vector(u.NA)) == TRUE)) > 
            0) {
            warning("Warning: NA in matrix U", call. = FALSE)
        }
    }
    matrix.u[u.NA] <- 0 
    matrix.X <- matrix.w %*% matrix.u
    return(list(matrix.w = matrix.w, matrix.u = matrix.u, matrix.X = matrix.X))
}