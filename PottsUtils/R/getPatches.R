getPatches <- function(bonds, nvertex){
    
    if(is.matrix(bonds) && ncol(bonds) == 2){
        a <- bonds[,1] - 1
        b <- bonds[,2] - 1
        nbond <- nrow(bonds)
    }
    else
        if(is.vector(bonds) && length(bonds) == 2){
            a <- bonds[1] - 1
            b <- bonds[2] - 1
            nbond <- 1
        }
        else stop("'bonds' has to be a matrix of 2 columns or a vector of length 2.")
    if(!(length(bonds)==0) && nvertex < max(bonds))
        stop("The total number of vertices is less than the largest vertex number in 'bonds'.")

    a <- as.integer(a)
    b <- as.integer(b)
    patches <- .Call("getPatches", a, b, nbond, nvertex)
    
    lapply(1:length(patches), function(i) patches[[i]]+1)
    
}

