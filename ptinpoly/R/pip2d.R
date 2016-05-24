pip2d = function(Vertices,Queries) {
    # Output is a NUMERIC vector.

    # Basic checks of Vertices input argument.
    vertDims = dim(Vertices)
    numColsV = vertDims[2]
    if ( numColsV != 2 ) {
        print("pip2d(): Number of columns in Vertices must be 2!")
        return(as.vector(rep(-3,nrow(Queries))))
    }

    # Basic checks of Queries input argument.
    querDims = dim(Queries)
    numColsQ = querDims[2]
    if ( numColsQ != 2 ) {
        print("pip2d(): Number of columns in Queries must be 2!")
        return(as.vector(rep(-6,nrow(Queries))))
    }

    # Invoke the ptinpoly C++ code (for convenience, it has
    # the same name as this R function, but this wasn't necessary)
    Output =.C("pip2d",
        vts    = as.double(Vertices),
        nVts   = as.integer(nrow(Vertices)),
        qrs    = as.double(Queries),
        nQrs   = as.integer(nrow(Queries)),
        Result = as.vector(rep(0,nrow(Queries)),mode="integer")
    )

    # Return the Results data.
    return(Output$Result)
}
