#' Calculate Dirichlet normal energy of a surface
#'
#' A function that calculates Dirichlet normal energy following the method of Bunn et
#' al. (2011) Comparing Dirichlet normal surface energy of tooth crowns, a new
#' technique of molar shape quantification for dietary inference, with previous methods
#' in isolation and in combination. Am J Phys Anthropol 145:247-261 doi: 10.1002
#' ajpa.21489
#'
#' @param plyFile An object of class 'mesh3d' and 'shape3d' with calculated normals
#' @param outliers The percentile of Dirichlet energy density values to be excluded 
#' defaults to top 0.1 percent
#' @param BoundaryDiscard Logical indicating how to handle the exclusion of 
#' boundary faces. Defaults to Leg which exlcudes faces which have a leg on the
#' boundary
#'
#' @details The function requires an object created by reading in a ply file utilizing
#' either the read.ply or the read.AVIZO.ply function, with calculated normals.
#'
#' Dirichlet normal energy is calculated on meshes that represent specimen surfaces and
#' have already been simplified to 10,000 faces and pre-smoothed in a 3D data
#' editing program. 
#'
#' In the default settings, the function seeks to discard boundary faces. This can be 
#' changed by adjusting the BoundaryDiscard argumen to 'None' which will not discard 
#' any faces on on the boundary. Further, there are two ways of excluding boundary faces.
#' Either if they have a leg on the boundary by setting BoundaryDiscard='Leg' or by 
#' excluding any face which has a vertex on the boundary with BoundaryDiscard='Vertex'.
#' The function defaults to remove the top 0.1 percent of calculated energy densities as 
#' outliers. Mesh orientation does not affect for this calculation.
#'
#' @importFrom
#' stats quantile aggregate
#'
#' @importFrom
#' Rvcg vcgGetEdge
#'
#' @export
#' DNE





DNE <- function(plyFile, outliers=0.1, BoundaryDiscard='Leg') {
	if(BoundaryDiscard!='Leg' && BoundaryDiscard!='Vertex' && BoundaryDiscard!='None'){
		stop("BoundaryDiscard must be set to 'Leg' 'Vertex' or 'None'. Selected 'None' if you are working with a closed surface.")
	}
	size <- cSize(plyFile$vb[-4,])
	plyFile$vb <- plyFile$vb/size*100
	
	ply <- Equal_Vertex_Normals(plyFile) ## Correct the Vertex Normals Calculation
	Es <- compute_energy_per_face(ply) ## Compute DNE values for each face of the surface
	
	if(BoundaryDiscard!='None'){
	### Extracting and removing Edge Facess
	if(BoundaryDiscard=='Leg') {
		edges <- vcgGetEdge(plyFile)
		temp <- subset(edges, edges$border==1)
		EdgeEs <- sort(as.numeric(temp$facept))
	}
	if(BoundaryDiscard=='Vertex') {
		edges <- vcgGetEdge(plyFile)
		bounds <- subset(edges, edges$border==1)
		edgeverts <- unique(c(bounds$vert1, bounds$vert2))
		list <- vertex_to_face_list(plyFile)
		EdgeEs <- sort(unique(unlist(list[edgeverts])))
	}
	
	Boundary_Values <- Es[EdgeEs,] ## This to be Exported, values of edge faces
	
	Es[EdgeEs,]$Dirichlet_Energy_Densities <- 0
	}
	
	### Extracting and removing outliers
	outs <- (100-outliers)/100
	
	DNEs <- Es$Dirichlet_Energy_Densities
	FAs <- Es$Face_Areas
	Q <- quantile(DNEs, probs=c(outs))
	
	Outlier_List <- which(DNEs > Q)
	
	Outliers <- Es[Outlier_List,] ## This to be Exported, outlier values
	
	DNEs[Outlier_List] <- 0
	
	CleanEs <- data.frame(Dirichlet_Energy_Densities=DNEs, Face_Areas=FAs)
	
	Surface_DNE <- sum(CleanEs$Dirichlet_Energy_Densities*CleanEs$Face_Areas)
	
	CleanEs[,1] <- CleanEs[,1]/size*100
	CleanEs[,2] <- CleanEs[,2]*size/100
	
	Outliers[,1] <- Outliers[,1]/size*100
	Outliers[,2] <- Outliers[,2]*size/100
	
	if(BoundaryDiscard!='None'){
	Boundary_Values[,1] <- Boundary_Values[,1]/size*100
	Boundary_Values[,2] <- Boundary_Values[,2]*size/100
	}
	
	nor <- ply$vb[4,]
	vbs <- plyFile$vb*size/100
	vbs <- vbs[-4,]
	rebuild <- rbind(vbs, nor)
	plyFile$vb <- rebuild
	
	if(BoundaryDiscard=='None') {
		Boundary_Values="No Boundary or No Boundary Discard Selected"
	}
	
	Out <- list(Surface_DNE=Surface_DNE, Face_Values=CleanEs, Boundary_Values=Boundary_Values, Outliers=Outliers, "plyFile"=plyFile)
	cat("Total Surface DNE =", Surface_DNE)
	return(Out)
	
}