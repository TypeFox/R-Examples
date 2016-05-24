"difference.vector" <-
  function(xyz, xyz.inds=NULL, normalize=FALSE) {

  xyz <- as.matrix(xyz)
  if (dim(xyz)[1L] < 2)
    stop("xyz must be a matrix with two rows")
  if (dim(xyz)[2L] < 6)
    stop("xyz does not contain sufficient coordinates")
  if (dim(xyz)[1L] > 2) {
    xyz <- xyz[1:2,]
    warning("xyz has more than two rows - using only the two first")
  }
  
  if ( is.null(xyz.inds) )
    xyz.inds <- seq(1, ncol(xyz))
  
  if ( length(which(is.na(xyz[,xyz.inds]))) > 0 )
    stop("xyz has NA values")
  
  a <- xyz[1, xyz.inds]
  b <- xyz[2, xyz.inds]

  if (length(a)!=length(b))
    stop("unequal lengths of the two coordinate sets")

  diff <- b-a
  if(normalize)
    diff <- normalize.vector(diff)
  
  return( diff )
}
