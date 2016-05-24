##' circulantij function
##'
##' A function to return the "idx" i.e. c(i,j) element of a circulant matrix with base "base".
##'
##' @param idx vector of length 2 th (i,j) (row,column) index to return 
##' @param base the base matrix of a circulant matrix
##' @return the ij element of the full circulant
##' @export
circulantij <- function(idx,base){
    i <- idx[1]
    j <- idx[2]
    m <- nrow(base)
    n <- ncol(base)
    if(i%%m==0){
        I1 <- i/m
    }
    else{        
        I1 <- floor(i/m) + 1
    }
    
    if(j%%m==0){
        J1 <- j/m
    }
    else{
        J1 <- floor(j/m) + 1
    }
    colno <- (J1-I1)%%n + 1
    
    I2 <- i%%m + 1
    J2 <- j%%m + 1
    rowno <- (J2-I2)%%m + 1

    return(base[rowno,colno])
}


##' circulant function
##'
##' generic function for constructing circulant matrices
##'
##' @param x an object
##' @param ... additional arguments
##' @return method circulant
##' @export
circulant <- function(x,...){
    UseMethod("circulant")
}



##' circulant.numeric function
##'
##' returns a circulant matrix with base x
##'
##' @method circulant numeric
##' @param x an numeric object
##' @param ... additional arguments
##' @return a circulant matrix with base x
##' @export
circulant.numeric <- function(x,...){
    n <- length(x)
    M <- matrix(NA,n,n)
    idx <- 1:n
    for(i in 1:n){
        M[i,] <- x[idx]
        idx <- c(idx[n],idx[-n])
    }
    return(M)
}



##' circulant.matrix function
##'
##' If x is a matrix whose columns are the bases of the sub-blocks of a block circulant matrix, then this function returns the 
##' block circulant matrix of interest.
##'
##' @method circulant matrix
##' @param x a matrix object
##' @param ... additional arguments
##' @return If x is a matrix whose columns are the bases of the sub-blocks of a block circulant matrix, then this function returns the block circulant matrix of interest.
##' @export
circulant.matrix <- function(x,...){
    M <- dim(x)[1]
    N <- dim(x)[2]
    submats <- t(apply(x,2,circulant))
    mat <- matrix(NA,M*N,M*N)
    idx <- 1:N
    ct <- 1
    for (i in 1:N){
        for (j in 1:N){ 
            xstart <- M*floor((ct-1)/N) + 1
            mult <- ct%%N
            if (mult==0){
                mult <- N
            }
            ystart <- M*(mult-1) + 1
            mat[xstart:(xstart+M-1),ystart:(ystart+M-1)] <- submats[idx[j],]
            ct <- ct + 1
        }
        idx <- c(idx[N],idx[-N])
    }
    return(mat)
}