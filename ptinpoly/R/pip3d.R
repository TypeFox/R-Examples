pip3d = function(Vertices,Faces,Queries) {
    # Output is a NUMERIC vector.

    # Basic checks of Vertices input argument.
    vertDims = dim(Vertices)
    numColsV = vertDims[2]
    if ( numColsV != 3 ) {
        print("pip3d(): Number of columns in Vertices must be 3!")
        return(as.vector(rep(-3,nrow(Queries))))
    }

    # Basic checks of Faces input argument.
    faceDims = dim(Faces)
    numColsF = faceDims[2]
    if ( numColsF != 3 ) {
        print("pip3d(): Number of columns in Faces must be 3!")
        return(as.vector(rep(-4,nrow(Queries))))
    }

    if ( min(Faces) < 1 ) {
        print("pip3d(): Values in Faces must be greater than 0!")
        return(as.vector(rep(-5,nrow(Queries))))
    }

    # Check to ensure that ALL vertices in the matrix VERTICES
    # are actually used in the matrix defining the faces, FACES.
    # First, make a list of the vertices in VERTICES that
    # are ACTUALLY USED in FACES.
    numVertices <- nrow(Vertices)
    listOfUsedVertices <- unique(as.numeric(Faces))
    if ( ! all(1:numVertices %in% listOfUsedVertices) ) {
        # If we got in this IF block, then not all vertices
        # in the matrix VERTICES are actually used in FACES.
        # Modify these two matrices to so that ALL vertices
        # listed in the matrix VERTICES are actually used
        # in the matrix FACES, and that the indices in FACES
        # run from exactly 1 to exactly the number of vertices
        # defined in the matrix VERTICES.

        # Make an updated version of Vertices that contains
        # ONLY the vertices that are mentioned in FACES.
        Vertices <- Vertices[listOfUsedVertices,]

        # Make an updated version of FACES that references
        # the vertices in the updated Vertices.
        numVerticesUsed <- length(listOfUsedVertices)
        FacesTmp <- matrix(NA,nrow(Faces),ncol(Faces)); # Initialize to NA.
        for ( i in 1:numVerticesUsed ) {
            vertIndex <- listOfUsedVertices[i]
            FacesTmp[ Faces == vertIndex ] <- i
        }
        Faces <- FacesTmp
    }

    # Basic checks of Queries input argument.
    querDims = dim(Queries)
    numColsQ = querDims[2]
    if ( numColsQ != 3 ) {
        print("pip3d(): Number of columns in Queries must be 3!")
        return(as.vector(rep(-6,nrow(Queries))))
    }

    # Invoke theptinpoly C++ code (for convenience, it has
    # the same name as this R function, but this wasn't necessary)
    Output =.C("pip3d",
        vts    = as.double(Vertices),
        nVts   = as.integer(nrow(Vertices)),
        fcs    = as.integer(Faces),
        nFcs   = as.integer(nrow(Faces)),
        qrs    = as.double(Queries),
        nQrs   = as.integer(nrow(Queries)),
        Result = as.vector(rep(0,nrow(Queries)),mode="integer")
    )

    # Return the Results data.
    return(Output$Result)
}
