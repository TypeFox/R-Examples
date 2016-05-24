#' Calculate Boyer's (2008) relief index for a surface
#'
#' A function that calculates relief index following Boyer (2008) Relief index of
#' second mandibular molars is a correlate of diet among prosimian primates and
#' other mammals. J Hum Evol 55:1118-1137 doi: 10.1016/j.jhevol.2008.08.002
#'
#' @param plyFile An object of classes 'mesh3d' and 'shape3d'
#' @param alpha Step size for calculating the outline. See details
#'
#' @details The function requires an object created by reading in a ply file utilizing
#' either the read.ply or the read.AVIZO.ply function, with calculated normals.
#' 
#' Relief index is calculated by the ratio of three-dimensional surface area to two
#' dimensional area on meshes that represent specimen surfaces and have already
#' been pre-smoothed in a 3D data editing program. Alignment of the point cloud
#' will have a large effect on patch orientation and must be done in a 3D data editing
#' program or {auto3dgm} prior to creating and reading in the ply file. The mesh
#' must be oriented such that the occlusal plane is parallel to the X- and Y-axes and
#' perpendicular to the Z-axis.
#'
#' Some files may fail with pancake[TempF,] : subscript out
#' of bounds. In these files it may be necessary to increase
#' the alpha value which is default set to 0.01. Increasing the
#' alpha value can cause the RFI function to over-estimate
#' the size of the footprint. Caution should be exercised when
#' troubleshooting by adjusting alpha
#'
#' @importFrom
#' alphahull ahull
#'
#' @export
#' RFI

RFI <- function(plyFile, alpha=0.01) {
  
  size <- cSize(plyFile$vb[-4,])
  plyFile$vb <- plyFile$vb/size*100
  
  ThreeDVerts1 <- t(plyFile$vb) ## Read in and properly transform original 3D vertices
	ThreeDVerts2 <- ThreeDVerts1[,-4]
	ThreeDFaces <- t(plyFile$it) ## Read in and properly transform the faces
	
	###### The Following section calculates the 3D surface area
	ThreeDFace_areas <- numeric(length(ThreeDFaces[,1])) ## Repository vector for 3D face areas
	
	for (i in 1:length(ThreeDFace_areas)) {
		TempF <- ThreeDFaces[i,] ## Pull each face one at a time
		TempV <- ThreeDVerts2[TempF,] ## Pull each vert from the designated face
		
		b1 <- TempV[2,]-TempV[1,] ## Begin Calculations
		b2 <- TempV[3,]-TempV[1,]
		g <- matrix(c(sum(b1*b1), sum(b1*b2), sum(b2*b1), sum(b2*b2)), nrow=2)
		
		ThreeDFace_areas[i] <- 0.5*sqrt(abs(g[1,1]*g[2,2]-g[1,2]*g[2,1]))
	}
	ThreeDArea <- sum(ThreeDFace_areas)*size/100 ## This will be an essential export 
	
	
	###### The Following sectionc calculates the 2D surface area
	x <- ThreeDVerts2[,1] - mean(ThreeDVerts2[,1])
	y <- ThreeDVerts2[,2] - mean(ThreeDVerts2[,2])
	z <- ThreeDVerts2[,3] - mean(ThreeDVerts2[,3])
	
	Shifted <- as.matrix(cbind(x, y, z)) ## Centered pt cloud at the origin. 
	Origin <- c(0,0,0)
	Shifted <- rbind(Shifted, Origin) ## Added the origin pt. 
	
	pancake <- as.matrix(cbind(Shifted[,1:2], z=rep(0, length(Shifted[,1]))))
	
	hull <- ahull(pancake[,1:2], alpha=alpha)  ## calculate alpha hull which rings the flattened point cloud
	
	arcs <- hull$arcs ## Begin building pie-slice triangles
	if(length(arcs[,7])!=length(unique(arcs[,7]))){
		stop('alpha adjustment required')
	}
	STedges <- arcs[,'end1']
	EDedges <- arcs[,'end2']
	test <- STedges-EDedges
	if(length(which(test==0))>0) {
		stop('alpha adjustment required')
	}
	Or <- length(Shifted[,1])
	center <- rep(Or, length(STedges))
	slices <- cbind(center, STedges, EDedges)
	
	TwoDFace_areas <- numeric(length(slices[,1])) ## Repository vector for 2D face areas
	for (i in 1:length(TwoDFace_areas)) {
		TempF <- slices[i,] ## Pull each slice one at a time
		TempV <- pancake[TempF,] ## Pull each flattened vert from the designated face
		
		b1 <- TempV[2,]-TempV[1,] ## Begin Calculations
		b2 <- TempV[3,]-TempV[1,]
		g <- matrix(c(sum(b1*b1), sum(b1*b2), sum(b2*b1), sum(b2*b2)), nrow=2)
		
		TwoDFace_areas[i] <- 0.5*sqrt(abs(g[1,1]*g[2,2]-g[1,2]*g[2,1]))
	}
	TwoDArea <- sum(TwoDFace_areas)*size/100 ## This will be an essential export
	
	RFI <- log(sqrt(ThreeDArea)/sqrt(TwoDArea))  
  Out <- list(Surface_RFI=RFI, Three_D_Area=ThreeDArea, Two_D_Area=TwoDArea, Translated_Pts=Shifted, Flattened_Pts=pancake, Footprint_Triangles=slices, "plyFile"=plyFile)
  cat("RFI =", RFI, "\n")
  cat("3D Area =", ThreeDArea, "\n")
  cat("2D Area =", TwoDArea)
  return(Out)
}