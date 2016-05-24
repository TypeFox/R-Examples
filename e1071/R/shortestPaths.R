allShortestPaths <- function(x){
    
    x <- as.matrix(x)
    x[is.na(x)] <- .Machine$double.xmax
    x[is.infinite(x) & x>0] <- .Machine$double.xmax
    if(ncol(x) != nrow(x))
        stop("x is not a square matrix")
    n <- ncol(x)
    
    z <- .C("e1071_floyd",
            as.integer(n),
            double(n^2),
            as.double(x),
            integer(n^2),
            PACKAGE = "e1071")
    z <- list(length = matrix(z[[2]], n),
              middlePoints = matrix(z[[4]]+1, n))
    z$length[z$length == .Machine$double.xmax] <- NA
    z
}
                       
extractPath <- function(obj, start, end){

    z <- integer(0)
    path <- function(i, j){

        k <- obj$middlePoints[i, j]
        if (k != 0)
        {
            path(i,k);
            z <<- c(z, k)
            path(k,j);
        }
    }
    path(start,end)
    c(start, z, end)
}
    
