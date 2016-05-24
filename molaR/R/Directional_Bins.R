#' This bins the faces into directional categories
#'
#' bins into 8 directional categories on the basis of their orientations
#' @param plyFile a stanford PLY file 
#' @param rotation the amount to rotate the specimen by
#' Directional_Bins()


Directional_Bins <- function(plyFile, rotation=0) {

FNormals <- plyFile$Face_Normals
rotation <- rotation*pi/180

Direction <- numeric(length(FNormals[1,]))

if(rotation==0) {
for (i in 1:length(Direction)) {
	
	pts <- FNormals[,i]
	Direction[i] <- atan2(pts[2], pts[1])
	}
}
if(rotation!=0) {
	for (i in 1:length(Direction)) {
		pts <- FNormals[,i]
		xRot <- pts[1]*cos(rotation) - pts[2]*sin(rotation)
		yRot <- pts[1]*sin(rotation) + pts[2]*cos(rotation)
		Direction[i] <- atan2(yRot, xRot)
	}
}

DBins <- .bincode(Direction, c(-pi, -3/4*pi, -1/2*pi, -1/4*pi, 0, 1/4*pi, 1/2*pi, 3/4*pi, pi))
bads <- which(is.na(DBins))
for(i in bads) {
	replace = i+1
	while(is.na(DBins[replace])){
		if(replace+1 > length(DBins)) replace=1
		replace=replace+1
	}
	DBins[i] <- DBins[replace]
}

plyFile$Directional_Bins <- DBins
return(plyFile)
}