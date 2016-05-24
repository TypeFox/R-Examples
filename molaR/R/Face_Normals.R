#' Function to find Face Normals
#'
#' This function re-computes the face normals in a way consistent with MorphoTester.
#' @param plyFile a stanford PLY file  
#' Face_Normals()

Face_Normals <- function(plyFile){

Faces <- plyFile$it

plyFile$Face_Normals <- matrix(0, nrow=length(Faces[1,]), ncol=3)


FNormals <- plyFile$Face_Normals

verts <- plyFile$vb


for (i in 1:length(FNormals[,1])) {
	
	pts <- Faces[,i]
	
	pt1 <- verts[,pts[1]]
	pt2 <- verts[,pts[2]]
	pt3 <- verts[,pts[3]]
	
	Vec1 <- pt2 - pt1
	Vec2 <- pt3 - pt1
	
	pfNorm <- as.vector(c(Vec1[2]*Vec2[3]-Vec1[3]*Vec2[2], Vec1[3]*Vec2[1]-	Vec1[1]*Vec2[3],Vec1[1]*Vec2[2]-Vec1[2]*Vec2[1]))
	Length <- sqrt(sum(pfNorm^2)) 
	FNormals[i,] <- pfNorm/Length
	
}

plyFile$Face_Normals <- t(FNormals)
return(plyFile)

}
