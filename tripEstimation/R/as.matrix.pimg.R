`as.matrix.pimg` <-
function(x, ...) {

  pimg <- x
  img <- matrix(0,pimg$xbound[3],pimg$ybound[3])
  if(!is.null(pimg$image)) {
    off <- pimg$offset
    img[off[1]:(off[1]+nrow(pimg$image)-1),
        off[2]:(off[2]+ncol(pimg$image)-1)] <- pimg$image
  }
  img
}

