#' Important function for re-doing the vertex normals for the DNE calculation. 
#'
#' The geomorph import function does not generate the correct vertex normals.
#' @param plyFile a stanford PLY file  
#' Equal_Vertex_Normals()


Equal_Vertex_Normals <- function(plyFile) {
	
	VertFace <- vertex_to_face_list(plyFile) ## list of Vertices with associated faces
	FaceVert <- plyFile$it ## List of Faces with assocated vertices, the order follow right hand rule. 
	v <- plyFile$vb ## list of vertices in space
	
	rawNorms <- plyFile$normals * 0  ## Set up empty matrix for storing summed normals for each vertex. 
	rownames(rawNorms) <- c('x', 'y', 'z', 'length')

	for (i in 1:length(VertFace)) {

	Faces <- VertFace[[i]]
		VertFaceMatrix <- matrix(0, nrow=length(Faces), ncol=3) ## Empty matrix for storing normal of each face at a vertex
		colnames(VertFaceMatrix) <- c('xpts', 'ypts', 'zpts')

		## The following for loop produces normals on the given vertex from each of the faces on that vertix. It assumes that the ply file face-to-vertex listing follows right hand rule. Such that points 1 2 3 follow counter clockwise. This is standard for STL- and ply files. 
		for (j in 1:length(Faces)) {
			Face <- Faces[j]
			pts <- FaceVert[,Face]

			pt1 <- v[,pts[1]]
			pt2 <- v[,pts[2]]
			pt3 <- v[,pts[3]]

				if (which(pts==i) == 1) {
					Vec1 <- pt2 - pt1
					Vec2 <- pt3 - pt1
						}
				if (which(pts==i) == 2) {
					Vec1 <- pt3 - pt2
					Vec2 <- pt1 - pt2
						}
				if (which(pts==i) == 3) {
					Vec1 <- pt1 - pt3
					Vec2 <- pt2 - pt3
						}
						
			## The following is the cross-product for two vectors. Following right-hand rule this cross product requires that Vec2 be counterclockwise to Vec1. STL and ply files are arranged as such. 
			pfNorm <- as.vector(c(Vec1[2]*Vec2[3]-Vec1[3]*Vec2[2], Vec1[3]*Vec2[1]-	Vec1[1]*Vec2[3],Vec1[1]*Vec2[2]-Vec1[2]*Vec2[1]))
			Length <- sqrt(sum(pfNorm^2)) 
			VertFaceMatrix[j,] <- pfNorm/Length
			}
		vNormRaw <- as.vector(colSums(VertFaceMatrix))
		vNormRaw <- as.vector(c(vNormRaw, sqrt(sum(vNormRaw^2))))
		rawNorms[,i] <- vNormRaw
	}

Norms <- rawNorms %*% diag(1/rawNorms[4,]) ## Normalize, this sets all to unit length

plyFile$normals <- Norms
return(plyFile)
}
