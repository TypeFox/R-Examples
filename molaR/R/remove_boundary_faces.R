#' Remove boundary faces
#' 
#' Important function for masking the edge faces
#' @param plyFile a stanford PLY file 
#' @param Energy_Per_Face_Values information on E per face
#' remove_boundary_faces()

remove_boundary_faces <- function(Energy_Per_Face_Values, plyFile) {
	
	EdgeVerts <- edge_vertices(plyFile) ## identify vertices on the edge
	Vert_to_Face <- vertex_to_face_list(plyFile) ## produce list of vertices and their associated faces
	
	Es_to_be_dropped <- sort(unique(unlist(Vert_to_Face[EdgeVerts])))
	
	Energy_Per_Face_Values[Es_to_be_dropped,]$DNE_Values <- 0
	
	return(Energy_Per_Face_Values)
}
