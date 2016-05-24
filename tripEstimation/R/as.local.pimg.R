"as.local.pimg" <-
function(pimg) {

  img <- coords.pimg(pimg)
  img$x <- img$x[pimg$offset[1]:(pimg$offset[1] + nrow(pimg$image)-1)]
  img$y <- img$y[pimg$offset[2]:(pimg$offset[2] + ncol(pimg$image)-1)]
  img$z <- pimg$image
  img
}

