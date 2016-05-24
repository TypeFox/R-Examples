permutations <- function(n){
   ## from package e1071
    if (n == 1) 
        return(matrix(1))
    else if (n < 2) 
        stop("n must be a positive integer")
    z <- matrix(1)
    for (i in 2:n) {
        x <- cbind(z, i)
        a <- c(1:i, 1:(i - 1))
        z <- matrix(0, ncol = ncol(x), nrow = i * nrow(x))
        z[1:nrow(x), ] <- x
        for (j in 2:i - 1) {
            z[j * nrow(x) + 1:nrow(x), ] <- x[, a[1:i + j]]
        }
    }
    dimnames(z) <- NULL
    z
}

invperm <- function(perm){
        sort(perm,index.return=TRUE)$ix
    }

reord <- function(hilf, perm){
    ## reorder matrix of binary numbers
    ## generating factors are isomorphic with binary numbers
    ## A=1, B=10, C=100, D=1000, usw.
    ##     AB=A+B=11 usw.

    ## column numbers for reordered generating factors can be obtained by 
    ##     switching digit positions

    ## reord does this and prepares the resulting matrix for calculating the column numbers
    ## with package sfsmisc
    aus <- hilf[nrow(hilf):1,,drop=FALSE][perm,,drop=FALSE][nrow(hilf):1,,drop=FALSE]
    class(aus) <- "basedInt"
    attr(aus,"base") <- 2
    aus
  }


getNext <- function(perm){
   ## function for next permutation in lexicographical order
   ## adapted from http://www.cut-the-knot.org/do_you_know/AllPerm.shtml
   ##     provided by Alexander Bogomolny, based on Dijkstra 1997 p.71
    N <- length(perm)
    i <- N

    while (perm[i-1] >= perm[i] & i>2) i <- i-1
 
    j <- N+1

    while (perm[j-1] <= perm[i-1] & j>2) {
         j <- j-1
       }

   ## swap values at positions (i-1) and (j-1)
    perm[c(i-1, j-1)] <- perm[c(j-1, i-1)]

    i <- i+1; j <- N+1

    while (i < j)
    {
      perm[c(i-1, j-1)] <- perm[c(j-1, i-1)]
      i <- i+1; j <- j-1
    }
    perm
}
