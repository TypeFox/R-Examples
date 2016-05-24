matrix.p<-function (comm, dist.spp, notification = TRUE)
{
	comm<-as.matrix(comm)
	dist.spp<-as.matrix(dist.spp)
    matrix.w <- sweep(comm, 1, rowSums(comm, na.rm=TRUE), "/")
	w.NA <- apply(matrix.w, 2, is.na)
    matrix.w[w.NA] <-0
    if(notification==TRUE){
    	if(length(which(unique(as.vector(w.NA))==TRUE))>0)
    	{
			warning("Warning: NA in community data",call.=FALSE)		
    	}  	 
    }
    similar.phy <- 1 - (dist.spp/max(dist.spp))
    matrix.phy <- 1/colSums(similar.phy)
    matrix.q <- sweep(similar.phy, 1, matrix.phy, "*")
    matrix.P <- matrix.w %*% matrix.q
    return(list(matrix.w = matrix.w, matrix.q = matrix.q, matrix.P = matrix.P))
}