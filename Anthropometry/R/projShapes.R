projShapes <- function(clust, array3D, asig, prototypes){
  out_proc <- c()
  out_proc <- shapes::procGPA(array3D[, , asig == clust], distances = TRUE, pcaoutput = TRUE)
  
  shapes::plotshapes(out_proc$rotated)
  points(prototypes[, , clust], col = 2)
}