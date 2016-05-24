#' Clean up problem ply files
#'
#' Function will remove floating verticies, and faces with zero
#' area. These can cause issues when using molaR's primary
#' functions of DNE, RFI, and OPC
#'
#' @param plyFile An object of classes 'mesh3d' and 'shape3d'
#' @param cleanType logical asking what to clean, Verticies,
#' Faces or Both. Defaults to Both.
#'
#' @details This function cleans up problematic ply files. Some
#' smoothed files will have faces of zero area, or floating
#' verticies. DNE and OPC cannot be calculated on these files.
#' Running the plys through this function will allow those
#' calculations to be made. 
#'
#'
#' @export
#' molaR_Clean


molaR_Clean <- function(plyFile, cleanType='Both'){
  if(cleanType != "Both" && cleanType != "Faces" && cleanType != "Vertices"){
    stop("cleanType must be set to either 'Faces', 'Vertices', or 'Both'.")
  }
  if(cleanType=='Both'){
    newPly <- face_areas(plyFile)
    Areas <- newPly$Face_Areas
    Zeroes <- which(Areas[]==0)
    cat("Removed", length(Zeroes), "faces with area = 0\n")
    if(length(Zeroes)>0){
      plyFile$it <- plyFile$it[,-Zeroes]
      cat("Indices of removed faces:", Zeroes, "\n")
    }
    Verts <- 1:ncol(plyFile$vb)
    VertList <- as.vector(plyFile$it)
    InFaces <- Verts %in% VertList
    NotIn <- which(InFaces==FALSE)
    cat("Removed", length(NotIn), "unreferenced vertices from mesh")
    if(length(NotIn)>0){
      for(i in 1:length(NotIn)){
        plyFile$vb <- plyFile$vb[,-NotIn[i]]
        plyFile$normals <- plyFile$normals[,-NotIn[i]]
        HighFaces <- plyFile$it > NotIn[i]
        plyFile$it[HighFaces] <- plyFile$it[HighFaces]-1
      }
      cat("\nIndices of removed vertices:", NotIn)
    }
  }
  if(cleanType=='Vertices'){
    Verts <- 1:ncol(plyFile$vb)
    VertList <- as.vector(plyFile$it)
    InFaces <- Verts %in% VertList
    NotIn <- which(InFaces==FALSE)
    cat("Removed", length(NotIn), "unreferenced vertices from mesh")
    if(length(NotIn)>0){
      for(i in 1:length(NotIn)){
        plyFile$vb <- plyFile$vb[,-NotIn[i]]
        plyFile$normals <- plyFile$normals[,-NotIn[i]]
        HighFaces <- plyFile$it > NotIn[i]
        plyFile$it[HighFaces] <- plyFile$it[HighFaces]-1
      }
      cat("\nIndices of removed vertices:", NotIn)
    }
  }
  if(cleanType=='Faces'){
    newPly <- face_areas(plyFile)
    Areas <- newPly$Face_Areas
    Zeroes <- which(Areas[]==0)
    cat("Removed", length(Zeroes), "faces with area = 0")
    if(length(Zeroes)>0){
      plyFile$it <- plyFile$it[,-Zeroes]
      cat("\nIndices of removed faces:", Zeroes)
    }
  }
  return(plyFile)
}