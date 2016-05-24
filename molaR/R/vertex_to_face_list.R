#' function for making a list of faces on each vertex
#'
#' crucial function for getting a list of faces which will gather the faces per vertex.
#' @param plyFile a stanford PLY file  
#' vertex_to_face_list()

vertex_to_face_list <- function(plyFile){
	Faces <- t(plyFile$it) ## input the Faces matrix, must be from the read.AVIZO.ply format. 
	
	vnum <- max(Faces)     ## Number of vertices
	fnum <- length(Faces[,1]) ## Number of faces
	
	out <- vector('list', vnum) ## empty list length of number of vertices
	
	
	#### Array the face names in the vertex rows
	for (i in 1:fnum) {
		n1 <- Faces[i,1] # indexing the vertices from each face
		n2 <- Faces[i,2]
		n3 <- Faces[i,3]
		
		out[[n1]] <- paste(out[[n1]], i, sep=",") # pasting in the face names to each vertex row
		out[[n2]] <- paste(out[[n2]], i, sep=",")
		out[[n3]] <- paste(out[[n3]], i, sep=",")
		
	}
	
	#### Remove the first comma from each vertex row
	for (i in 1:vnum) {
		out[[i]] <- substring(out[[i]], 2)
		
	}
	
	#### Convert to numbers in each row element
	for (i in 1:vnum) {
		out[[i]] <- as.numeric(unlist(strsplit(out[[i]], ",")))
	}
	Vs <- plyFile$vb
	names(out) <- seq(1:length(Vs[1,]))
	
	return(out)
}
