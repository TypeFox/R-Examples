"trim.xyz" <- function(xyz, row.inds=NULL, col.inds=NULL, ...) {
  xyz <- as.xyz(xyz)

  if(is.select(row.inds))
    row.inds <- row.inds$xyz
  if(is.select(col.inds))
    col.inds <- col.inds$xyz
    
  if(!is.null(row.inds))
    xyz <- xyz[row.inds, , drop=FALSE]
  if(!is.null(col.inds))
    xyz <- xyz[, col.inds, drop=FALSE]

  return(as.xyz(xyz))
}
