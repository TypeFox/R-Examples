#' Function will compute the DNE per face.
#'
#' This will generate each Dirichlet's normal energy for each triangular face on the surface.
#' @param plyFile a stanford PLY file 
#' compute_energy_per_face()
#'
#' 


compute_energy_per_face <- function(plyFile) {
	
	Vertecies <- t(plyFile$vb) ## Read in and properly transform vertices
	Verts <- Vertecies[,-4]  
	
	Normals <- t(plyFile$normals) ## Read in and properly transform normals
	Norms <- Normals[,-4]
	
	Faces <- t(plyFile$it)  ## Read in and properly transform faces
	
	e <- numeric(length(Faces[,1])) ## Repository column for DNE-values (one per face)
	
	face_area <- numeric(length(Faces[,1]))  ## Repository column for face areas (one per face)
	
	for (i in 1:length(e)) {
		
		TempF <- Faces[i,] ## Properly pull the faces, normals and vertices for each iteration
		TempV <- Verts[TempF,]
		TempN <- Norms[TempF,]
		
		b1 <- TempV[2,]-TempV[1,] ## Begin Calculations
		b2 <- TempV[3,]-TempV[1,]
		
		g <- matrix(c(sum(b1*b1), sum(b1*b2), sum(b2*b1), sum(b2*b2)), nrow=2)
		
		if(1/rcond(g,'1') > 10^5) { ## Subset out unsolvable matrices
			
			e[i] <- 0
			face_area[i] <- 1
		}
		
		else {
			c1 <- TempN[2,] - TempN[1,]
			c2 <- TempN[3,] - TempN[1,]
			
			grumple <- matrix(c(sum(c1*c1), sum(c1*c2), sum(c2*c1), sum(c2*c2)), nrow=2)
			
			e[i] <- tr(solve(g)%*%grumple)
			face_area[i] <- 0.5* sqrt(g[1,1]*g[2,2]-g[1,2]*g[2,1])
		}
		
	}
	
	Values_Per_Face <- data.frame(Dirichlet_Energy_Densities=e, Face_Areas=face_area)
	return(Values_Per_Face)
	
}