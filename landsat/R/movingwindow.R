movingwindow <-
function(x, kernel)
{
    # apply kernel as a moving window to x
    results <- x
    x <- as.matrix(x)

    mwoffset <- (nrow(kernel)-1)/2

    newmat <- matrix(NA, nrow=nrow(x), ncol=ncol(x))

    for(i in (1+mwoffset):(nrow(x)-mwoffset)) {
        for(j in (1+mwoffset):(ncol(x)-mwoffset)) {
            newmat[i, j] <- sum(kernel * x[(i-mwoffset):(i+mwoffset), (j-mwoffset):(j+mwoffset)])
        }
    }

    # return the same structure as the input values
    if(class(results) == "SpatialGridDataFrame")
        results@data[,1] <- as.vector(newmat)
    else if(is.data.frame(results))
        results <- data.frame(matrix(newmat, nrow=nrow(results), ncol=ncol(results)))
    else if(is.matrix(results))
        results <- matrix(newmat, nrow=nrow(results), ncol=ncol(results))

    results
}

