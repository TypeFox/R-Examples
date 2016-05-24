#' Function to calculate face areas. 
#'
#' This function calculates the area of each face on a ply file
#' @param plyFile a stanford PLY file 
#' face_areas()

face_areas <- function(plyFile){
	
	ThreeDVerts1 <- t(plyFile$vb) ## Read in and properly transform original 3D vertices
	ThreeDVerts2 <- ThreeDVerts1[,-4]
	ThreeDFaces <- t(plyFile$it) 

ThreeDFace_areas <- numeric(length(ThreeDFaces[,1])) ## Repository vector for 3D face areas
	
	for (i in 1:length(ThreeDFace_areas)) {
		TempF <- ThreeDFaces[i,] ## Pull each face one at a time
		TempV <- ThreeDVerts2[TempF,] ## Pull each vert from the designated face
		
		b1 <- TempV[2,]-TempV[1,] ## Begin Calculations
		b2 <- TempV[3,]-TempV[1,]
		g <- matrix(c(sum(b1*b1), sum(b1*b2), sum(b2*b1), sum(b2*b2)), nrow=2)
		
		ThreeDFace_areas[i] <- 0.5*sqrt(abs(g[1,1]*g[2,2]-g[1,2]*g[2,1]))
	}
plyFile$Face_Areas <- ThreeDFace_areas

return(plyFile)
}	

