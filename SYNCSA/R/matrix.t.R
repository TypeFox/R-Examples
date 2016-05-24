matrix.t<-function (comm, traits, scale = TRUE, notification = TRUE) 
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
    b.NA <- apply(traits, 2, is.na)
    traits[b.NA] <-0
    if(notification==TRUE){
    	if(length(which(unique(as.vector(b.NA))==TRUE))>0)
    	{
			warning("Warning: NA in traits matrix",call.=FALSE)	
    	}
    }
	matrix.b <- traits
    matrix.T <- matrix.w %*% traits
    if (scale == "TRUE") {
        matrix.traits <- apply(matrix.T^2, 2, sum)
        matrix.T <- sweep(matrix.T, 2, sqrt(matrix.traits), "/")
    }
    return(list(matrix.w = matrix.w, matrix.b = matrix.b, matrix.T = matrix.T))
}