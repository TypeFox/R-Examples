
## =============================================================================
##  Function to create isosurface triangles of an array
## =============================================================================

createisosurf <- function(x, y, z, colvar, level = mean(colvar, na.rm = TRUE))  {
  DD <- dim(colvar)
  if (length(x) != DD[1]) 
    stop ("'x' should be of length equal to first dimension of 'colvar'")
  if (length(y) != DD[2]) 
    stop ("'y' should be of length equal to second dimension of 'colvar'")
  if (length(z) != DD[3]) 
    stop ("'z' should be of length equal to third dimension of 'colvar'")

  Tri <- computeContour3d(vol = colvar, maxvol = max(colvar, na.rm = TRUE), 
     level = level, x = x, y = y, z = z, mask = NULL)
  colnames(Tri) <- c("x", "y", "z")
  invisible(Tri)
}  
