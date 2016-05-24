vegclustIndex <-function (y) {
	
    partition.coefficient <- function(clres) {
    	   if(inherits(clres,"vegclust")) clres = clres$memb
        xrows <- dim(clres)[1]
        partitioncoefficient <- sum(apply(clres^2, 1, sum))/xrows
        return(partitioncoefficient)
    }
    partition.entropy <- function(clres) {
    	   if(inherits(clres,"vegclust")) clres = clres$memb
        xrows <- dim(clres)[1]
        ncenters <- dim(clres)[2]
        partitionentropy <- 0
        for (i in 1:xrows) {
            for (k in 1:ncenters) {
                if (clres[i, k] != 0) 
                  partitionentropy <- partitionentropy + (clres[i,k] * log(clres[i, k]))
            }
        }
        partitionentropy <- partitionentropy/((-1) * xrows)
        return(partitionentropy)
    }
    if(inherits(y,"vegclust")) y = y$memb

    xrows <- dim(y)[1]
    ncenters <- dim(y)[2]
    findex = numeric(4)
    findex[1] = partition.coefficient(y)
    findex[2] = (findex[1]*ncenters-1)/(ncenters-1)
    findex[3] = partition.entropy(y)
    findex[4] = findex[3]/(1-(ncenters/xrows))
    names(findex)<-c("PC", "PCN","PE","PEN")
    return(findex)
}


