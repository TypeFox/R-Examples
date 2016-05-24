##
## imagematrix class definition
##
## $Header: /database/repository/rimage/R/Attic/imagematrix.R,v 1.1.2.5 2004/03/17 06:35:18 tomo Exp $
##
## Copyright (c) 2003 Nikon Systems Inc.
## For complete license terms see file LICENSE

imagematrix <- function(mat, type=NULL, ncol=dim(mat)[1], nrow=dim(mat)[2],
                        noclipping=FALSE) {
  if (is.null(dim(mat)) && is.null(type)) stop("Type should be specified.")
  if (length(dim(mat)) == 2 && is.null(type)) type <- "grey"
  if (length(dim(mat)) == 3 && is.null(type)) type <- "rgb"
  if (type != "rgb" && type != "grey") stop("Type is incorrect.")
  if (is.null(ncol) || is.null(nrow)) stop("Dimension is uncertain.")
  imgdim <- c(ncol, nrow, if (type == "rgb") 3 else NULL)
  if (length(imgdim) == 3 && type == "grey") {
    # force to convert grey image
    mat <- rgb2grey(mat)
  }
  if (noclipping == FALSE && ((min(mat) < 0) || (1 < max(mat)))) {
    warning("Pixel values were automatically clipped because of range over.") 
    mat <- clipping(mat)
  }
  mat <- array(mat, dim=imgdim)
  attr(mat, "type") <- type
  class(mat) <- c("imagematrix", class(mat))
  mat
}

print.imagematrix <- function(x, ...) {
  x.dim <- dim(x)
  cat("size: ", x.dim[1], "x", x.dim[2], "\n")
  cat("type: ", attr(x, "type"), "\n")
}

plot.imagematrix <- function(x, ...) {
  colvec <- switch(attr(x, "type"),
                grey=grey(x),
                rgb=rgb(x[,,1], x[,,2], x[,,3]))
  if (is.null(colvec)) stop("image matrix is broken.")
  colors <- unique(colvec)
  colmat <- array(match(colvec, colors), dim=dim(x)[1:2])
  image(x = 0:(dim(colmat)[2]), y=0:(dim(colmat)[1]),
        z = t(colmat[nrow(colmat):1, ]), col=colors,
        xlab="", ylab="", axes=FALSE, asp=1, ...)
}

imageType <- function(x) {
  attr(x, "type")
}

rgb2grey <- function(img, coefs=c(0.30, 0.59, 0.11)) {
  if (is.null(dim(img))) stop("image matrix isn't correct.")
  if (length(dim(img))<3) stop("image matrix isn't rgb image.")
  imagematrix(coefs[1] * img[,,1] + coefs[2] * img[,,2] + coefs[3] * img[,,3],
              type="grey")
}
clipping <- function(img, low=0, high=1) {
  img[img < low] <- low
  img[img > high] <- high
  img
}
normalize <- function(img) {
  (img - min(img))/(max(img) - min(img))
}

# the end of file
