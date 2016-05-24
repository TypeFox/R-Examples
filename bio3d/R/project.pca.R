"project.pca" <-
function(data, pca, angular=FALSE, fit=FALSE, ...) {
  
  if(angular)
    data <- wrap.tor(data)
  if(is.null(dim(data))) {
    if(ncol(pca$U) != length(data))
      stop("Dimensionality mismatch:  length(data)!=ncol(pca$U)")
    if(fit) data <- fit.xyz(pca$mean, data, ...)
    z <- (data - pca$mean) %*% pca$U
  } else {
    if(ncol(pca$U) != ncol(data))
      stop("Dimensionality mismatch:  ncol(data)!=ncol(pca$U)")
    if(fit) data <- fit.xyz(pca$mean, data, ...)
    z <- sweep(data, 2, pca$mean) %*% pca$U
  }
  return(z)
}

z2xyz.pca <- function(z.coord, pca) {

  
  if(is.null(nrow(z.coord))) {
    if( length(z.coord) > ncol(pca$U) )
      stop("Dimension miss-match: length(z.coord) > ncol(pca$U)")
    
    xyz.coord <- (z.coord  %*% t(pca$U[, c(1:length(z.coord)) ]) ) + pca$mean

  } else {
    if( ncol(z.coord) > ncol(pca$U) )
      stop("Dimension miss-match: ncol(z.coord) > ncol(pca$U)")

    xyz.coord <- NULL
    for(i in 1:nrow(z.coord)) {
      xyz.coord <- rbind(xyz.coord,
                         (z.coord[i,]  %*% t(pca$U[, c(1:length(z.coord[i,])) ]) ) + pca$mean)

    }
  }
  return(xyz.coord)
}


xyz2z.pca <- function(xyz.coord, pca) {
  return(project.pca(xyz.coord, pca))
}

