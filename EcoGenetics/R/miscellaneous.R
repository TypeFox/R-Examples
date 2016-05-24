#' Conversion of a non symmetric binary matrix into symmetric.
#' @param mat Matrix.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal


int.2symmetric <- function(mat) {
  mat[row(mat) > col(mat)] <- mat[row(mat) > col(mat)] + mat[row(mat) < col(mat)]
  mat[row(mat) > col(mat)][mat[row(mat) > col(mat)] > 0] <- 1
  mat[row(mat) < col(mat)] <- mat[row(mat) > col(mat)]
  mat
}


#' Creates a matrix without diagonal, in row order
#' @param mat Matrix.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

int.undimmattg <- function(mat) {
  ncolp <- ncol(mat) -1
  mat2 <- as.vector(t(mat))
  mat2 <-mat2[-which(col(mat) == row(mat))]
  mat2<-matrix(mat2, ncol = ncolp, byrow = T)
  mat2
}


#' Computing a distance matrix in meters among points in decimal degrees
#' under a spherical Earth model
#' @param XY data frame or matrix with latitude-longitude coordinates 
#' in decimal degrees format.
#' This program computes a distance matrix for Earth points in decimal degrees.
#' It assumes a spherical model with an Earth radius of 6371 km. 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal


int.dlatlon2distm <- function(XY) {
  out <- matrix(,nrow(XY), nrow(XY)) 
  for(i in 1:nrow(XY)) {
    for(j in 1:nrow(XY)) {
      lat1 <- XY[i, 1] 
      lon1 <- XY[i, 2] 
      lat2 <- XY[j, 1] 
      lon2 <- XY[j, 2] 
      R <- 6371                                
      dLat <- (lat2 - lat1) * pi / 180
      dLon <- (lon2 - lon1) * pi / 180
      a <- sin((dLat/2)) ^ 2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * (sin(dLon / 2)) ^2
      c <- 2 * atan2(sqrt(a), sqrt(1-a))
      d <- R * c  
      out[i, j] <- d
    }
  }
  rownames(out) <- rownames(XY)
  colnames(out) <- rownames(XY)
  as.dist(out)
}

