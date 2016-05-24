`as.image.pimg` <-
function(pimg) {
  img <- coords.pimg(pimg)
  img$z <- as.matrix.pimg(pimg)
  img
}

