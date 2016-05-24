########## R-function: bdiag ##########

# Builds from the list of the matrices 
# a block diagonal matrix

# Last changed 16 JUN 2006

bdiag <- function(x){
     if(!is.list(x)) stop("x not a list")
     n <- length(x)
     if(n==0) return(NULL)
     x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
stop("Zero-length component in x"))
     d <- array(unlist(lapply(x, dim)), c(2, n))
     rr <- d[1,]
     cc <- d[2,]
     rsum <- sum(rr)
     csum <- sum(cc)
     out <- array(0, c(rsum, csum))
     ind <- array(0, c(4, n))
     rcum <- cumsum(rr)
     ccum <- cumsum(cc)
     ind[1,-1] <- rcum[-n]
     ind[2,] <- rcum
     ind[3,-1] <- ccum[-n]
     ind[4,] <- ccum
     imat <- array(1:(rsum * csum), c(rsum, csum))
     iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
(y[3]+1):y[4]], imat=imat)
     iuse <- as.vector(unlist(iuse))
     out[iuse] <- unlist(x)
     return(out)
}
