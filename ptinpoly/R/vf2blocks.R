vf2blocks = function(Vertices,Faces) {
    # Conversion from Vertices-Faces representation to "block" representation:
    # first block for the coordinates of the first vertex of
    # each face, second block for the coordinates of the second vertex of
    # each face, and the third block for the coordinates of the third vertex
    # of each face.

    # Basic checks of Vertices input argument.
    vertDims = dim(Vertices)
    numColsV = vertDims[2] ;
	if ( numColsV != 3 ) {
	    print("vf2blocks(): Number of columns in Vertices must be 3!")
	    return(0);
	}

    # Basic checks of Faces input argument.
    faceDims = dim(Faces)
    numFaces = faceDims[1] ;
    numColsF = faceDims[2] ;
	if ( numColsF != 3 ) {
	    print("vf2blocks(): Number of columns in Faces must be 3!")
	    return(0);
	}
	if ( min(Faces) < 1 ) {
	    print("vf2blocks(): Values in Faces must be greater than 0!")
	    return(0);
	}

    # Loop over faces, build block matrices
    block1   = matrix(rep(0,3*numFaces),ncol = 3)
    block2   = matrix(rep(0,3*numFaces),ncol = 3)
    block3   = matrix(rep(0,3*numFaces),ncol = 3)
    for ( i in 1:numFaces ) {
        block1[i,] = Vertices[Faces[i,1],]
        block2[i,] = Vertices[Faces[i,2],]
        block3[i,] = Vertices[Faces[i,3],]
    }
    return(list(block1,block2,block3))
}
