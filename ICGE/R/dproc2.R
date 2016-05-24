dproc2 <- function(x, timepoints=NULL){
##### Distance by means of Procrustes statistic (two configurations) ##############################
# Input: 
# x: data matrix with non standarized profiles in rows 
# timepoints: considered timepoints. If timepoints=NULL, then timepoints=1:T
# Output: 
# d: procrustes distance
#############################################################################  
	x <- as.matrix(x)
	n <- dim(x)[1]
	p <- dim(x)[2]
	d2 <- matrix(0, nrow=n, ncol=n)
	for (i in 1:(n-1)){
		for (j in (i+1):n){
                        d2[i,j] <- proc2(x[i,], x[j,], timepoints=timepoints)   
			d2[j,i] <- d2[i,j]
		}
	}
        # Are there technical 0s?
        aux <- sum(d2<0)/2
        if (aux>0){
           cat(aux, " technical negative values are set to 0.")
           d2[d2<0] <- 0
        }
        
        d <- sqrt(d2)
        d <- as.dist(d)
	return(d)


} #function
