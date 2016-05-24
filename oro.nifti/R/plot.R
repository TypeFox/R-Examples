##
##
## Copyright (c) 2009-2011, Brandon Whitcher and Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
## $Id: plot.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## image() for class="nifti"
## promptMethods(image, "image-methods.Rd")
#############################################################################

image.nifti <- function(x, z=1, w=1, col=gray(0:64/64),
                        plane=c("axial", "coronal", "sagittal"),
                        plot.type=c("multiple", "single"), zlim=NULL,
                        xlab="", ylab="", axes=FALSE, oma=rep(0,4),
                        mar=rep(0,4), bg="black", ...) {
  switch(plane[1],
         "axial" = {
           aspect <- x@pixdim[3] / x@pixdim[2]
         },
         "coronal" = {
           if (length(dim(x)) == 3) {
             x@.Data <- aperm(x, c(1,3,2))
           } else {
             x@.Data <- aperm(x, c(1,3,2,4))
           }
           aspect <- x@pixdim[4] / x@pixdim[2]
         },
         "sagittal" = {
           if (length(dim(x)) == 3) {
             x@.Data <- aperm(x, c(2,3,1))
           } else {
             x@.Data <- aperm(x, c(2,3,1,4))
           }
           aspect <- x@pixdim[4] / x@pixdim[3]
         },
         stop(paste("Orthogonal plane", plane[1], "is not valid.")))
  ## set dimensions
  X <- nrow(x)
  Y <- ncol(x)
  Z <- nsli(x)
  W <- ntim(x)
  ## check dimensions
  if (is.na(X) || is.na(Y)) {
    stop("Missing rows/columns in NIfTI volume.")
  }
  if (! is.na(Z)) {
    if (z < 1 || z > Z) {
      stop("slice \"z\" out of range")
    }
  } else {
    plot.type <- "single"
  }
  ## check for z-limits in x; use internal by default
  if (is.null(zlim)) {
    zlim <- c(x@"cal_min", x@"cal_max")
    if (any(!is.finite(zlim)) || diff(zlim) == 0) {
      zlim <- c(x@"glmin", x@"glmax")
    }
    if (any(!is.finite(zlim)) || diff(zlim) == 0) {
      zlim <- range(x, na.rm=TRUE)
    }
  }
  breaks <- c(zlim[1], # min(x, zlim, na.rm=TRUE),
              seq(min(zlim, na.rm=TRUE), max(zlim, na.rm=TRUE),
                  length=length(col)-1),
              zlim[2]) # max(x, zlim, na.rm=TRUE))
  ## single or multiple images?
  if (plot.type[1] == "multiple") {
    index <- 1:Z
  } else {
    index <- z
  }
  lz <- length(index)
  ## plotting
  oldpar <- par(no.readonly=TRUE)
  par(mfrow=ceiling(rep(sqrt(lz),2)), oma=oma, mar=mar, bg=bg)
  if (is.na(Z)) { # two-dimensional matrix
    graphics::image(1:X, 1:Y, x, col=col, breaks=breaks, asp=aspect,
                    axes=axes, xlab=xlab, ylab=ylab, ...)
  } else {
    if (is.na(W)) { # three-dimensional array
      for (z in index) {
        graphics::image(1:X, 1:Y, x[,,z], col=col, breaks=breaks,
                        asp=aspect, axes=axes, xlab=xlab, ylab=ylab, ...)
      }
    } else { # four-dimensional array
      if (w < 1 || w > W)
        stop("volume \"w\" out of range")
      for (z in index) {
        graphics::image(1:X, 1:Y, x[,,z,w], col=col, breaks=breaks,
                        asp=aspect, axes=axes, xlab=xlab, ylab=ylab, ...)
      }
    }
  }
  par(oldpar)
  invisible()
}
#' @title Methods for Function `image'
#' 
#' @description Produce \dQuote{lightbox} layout of images for \code{nifti}, \code{anlz} and
#' \code{afni} objects.
#' 
#' @details Uses the S3 generic function \code{image}, with medical-image friendly
#' settings, to display \code{nifti}, \code{anlz} and \code{afni} class objects
#' in a \dQuote{lightbox} layout.
#' 
#' @name image-methods
#' @aliases image.nifti image-methods image,ANY-method image,afni-method
#' image,anlz-method image,nifti-method
#' @docType methods
#' @param x is an object of class \code{nifti} or similar.
#' @param z is the slice to be displayed (ignored when \code{plot.type =
#' "multiple"}).
#' @param w is the time point to be displayed (4D arrays only).
#' @param col is grayscale (by default).
#' @param plane is the plane of acquisition to be displayed (choices are
#' \sQuote{axial}, \sQuote{coronal}, \sQuote{sagittal}).
#' @param plot.type allows the choice between all slices being displayed, in a
#' matrix (left-to-right, top-to-bottom), or a single slice.
#' @param zlim is set to \code{NULL} by default and utilizes the internal image
#' range.
#' @param xlab is set to \dQuote{} since all margins are set to zero.
#' @param ylab is set to \dQuote{} since all margins are set to zero.
#' @param axes is set to \code{FALSE} since all margins are set to zero.
#' @param oma is the size of the outer margins in the \code{par} function.
#' @param mar is the number of lines of margin in the \code{par} function.
#' @param bg is the background color in the \code{par} function.
#' @param \dots other arguments to the \code{image} function may be provided
#' here.
#' @section Methods: \describe{ \item{x = "ANY"}{Generic function: see
#' \code{\link[graphics]{image}}.} \item{x = "nifti"}{Produce images for
#' \code{x}.} \item{x = "anlz"}{Produce images for \code{x}.} \item{x =
#' "afni"}{Produce images for \code{x}.} }
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{orthographic-methods}}, \code{\link{overlay-methods}}
#' @keywords methods
#' @import graphics
#' @import grDevices
#' @export
#' @rdname image-methods
setMethod("image", signature(x="nifti"), image.nifti)
#' @export
#' @rdname image-methods
setMethod("image", signature(x="anlz"), image.nifti)
#' @export
#' @rdname image-methods
setMethod("image", signature(x="afni"),
          function(x, ...) {
            x <- as(x, "nifti")
            image.nifti(x, ...)
          })

#############################################################################
## overlay() for class="nifti"
#############################################################################
#' Methods for Function overlay
#' 
#' Methods for function \code{overlay}
#' 
#' The \code{image} command is used multiple times to simultaneously visualize
#' one of the three orthogonal planes in two multidimensional arrays, one on
#' top of the other, for medical imaging data.
#' 
#' @name overlay-methods
#' @aliases overlay overlay-methods overlay,nifti,nifti-method
#' overlay,nifti,anlz-method overlay,nifti,array-method
#' overlay,anlz,nifti-method overlay,anlz,anlz-method overlay,anlz,array-method
#' overlay,array,nifti-method overlay,array,anlz-method
#' overlay,array,array-method overlay,afni,afni-method
#' overlay,afni,array-method overlay.nifti
#' @docType methods
#' @param x,y is an object of class \code{nifti} or similar.
#' @param z is the slice to be displayed (ignored when \code{plot.type =
#' "multiple"}).
#' @param w is the time point to be displayed (4D arrays only).
#' @param col.x is grayscale (by default).
#' @param col.y is hotmetal (by default).
#' @param zlim.x,zlim.y are set to \code{NULL} (by default) and taken from the
#' header information.
#' @param plane is the plane of acquisition to be displayed (choices are
#' \sQuote{axial}, \sQuote{coronal}, \sQuote{sagittal}).
#' @param plot.type allows the choice between all slices being displayed, in a
#' matrix (left-to-right, top-to-bottom), or a single slice.
#' @param xlab is set to \dQuote{} since all margins are set to zero.
#' @param ylab is set to \dQuote{} since all margins are set to zero.
#' @param axes is set to \code{FALSE} since all margins are set to zero.
#' @param oma is the size of the outer margins in the \code{par} function.
#' @param mar is the number of lines of margin in the \code{par} function.
#' @param bg is the background color in the \code{par} function.
#' @param \dots other arguments to the \code{image} function may be provided
#' here.
#' @section Methods: 
#' \describe{ 
#' \item{x = "nifti", y = "nifti"}{Produce overlay of \code{y} on \code{x}.} 
#' \item{x = "anlz", y = "anlz"}{Produce overlay of \code{y} on \code{x}.} 
#' \item{x = "afni", y = "afni"}{Produce overlay of \code{y} on \code{x}.} 
#' }
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{image-methods}}, \code{\link{overlay-methods}}
#' @keywords methods
#' @export
#' @rdname overlay-methods
overlay.nifti <- function(x, y, z=1, w=1, col.x=gray(0:64/64),
                          col.y=hotmetal(), zlim.x=NULL, zlim.y=NULL,
                          plane=c("axial", "coronal", "sagittal"),
                          plot.type=c("multiple","single"),
                          xlab="", ylab="", axes=FALSE, oma=rep(0,4),
                          mar=rep(0,4), bg="black", ...) {
  switch(plane[1],
         "axial" = {
           aspect <- x@pixdim[3] / x@pixdim[2]
         },
         "coronal" = {
           if (length(dim(x)) == 3) {
             x@.Data <- aperm(x, c(1,3,2))
           } else {
             x@.Data <- aperm(x, c(1,3,2,4))
           }
           y@.Data <- aperm(y, c(1,3,2))
           aspect <- x@pixdim[4] / x@pixdim[2]
         },
         "sagittal" = {
           if (length(dim(x)) == 3) {
             x@.Data <- aperm(x, c(2,3,1))
           } else {
             x@.Data <- aperm(x, c(2,3,1,4))
           }
           y@.Data <- aperm(y, c(2,3,1))
           aspect <- x@pixdim[4] / x@pixdim[3]
         },
         stop(paste("Orthogonal plane", plane[1], "is not valid.")))
  ## both volumes must have the same dimension
  if (! all(dim(x)[1:3] == dim(y)[1:3])) {
    stop("dimensions of \"x\" and \"y\" must be equal")
  }
  ## set dimensions
  X <- nrow(x)
  Y <- ncol(x)
  Z <- nsli(x)
  W <- ntim(x)
  ## check dimensions
  if (X == 0 || Y == 0 || Z == 0) {
    stop("size of NIfTI volume is zero, nothing to plot")
  }
  ## check for z-limits in x; use internal by default
  if (is.null(zlim.x)) {
    zlim.x <- c(x@"cal_min", x@"cal_max")
    if (any(!is.finite(zlim.x)) || diff(zlim.x) == 0) {
      zlim.x <- c(x@"glmin", x@"glmax")
    }
    if (any(!is.finite(zlim.x)) || diff(zlim.x) == 0) {
      zlim.x <- range(x, na.rm=TRUE)
    }
  }
  breaks.x <- c(min(x, zlim.x, na.rm=TRUE),
                seq(min(zlim.x, na.rm=TRUE), max(zlim.x, na.rm=TRUE),
                    length=length(col.x)-1),
                max(x, zlim.x, na.rm=TRUE))
  ## check for z-limits in y; use internal by default
  if (is.null(zlim.y)) {
    zlim.y <- c(y@"cal_min", y@"cal_max")
    if (any(!is.finite(zlim.y)) || diff(zlim.y) == 0) {
      zlim.y <- c(y@"glmin", y@"glmax")
    }
    if (any(!is.finite(zlim.y)) || diff(zlim.y) == 0) {
      zlim.y <- range(y, na.rm=TRUE)
    }
  }
  if (plot.type[1] == "multiple") {
    index <- 1:Z
  } else {
    index <- z
  }
  lz <- length(index)
  if (z < 1 || z > Z) {
    stop("slice \"z\" out of range")
  }
  oldpar <- par(no.readonly=TRUE)
  par(mfrow=ceiling(rep(sqrt(lz),2)), oma=oma, mar=mar, bg=bg)
  if (is.na(W)) { # three-dimensional array
    for (z in index) {
      graphics::image(1:X, 1:Y, x[,,z], col=col.x, breaks=breaks.x,
                      zlim=zlim.x, asp=aspect, axes=axes, xlab=xlab,
                      ylab=ylab, ...)
      graphics::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y, add=TRUE)
    }
  } else { # four-dimensional array
    if (w < 1 || w > W) {
      stop("volume \"w\" out of range")
    }
    for (z in index) {
      graphics::image(1:X, 1:Y, x[,,z,w], col=col.x, breaks=breaks.x,
                      zlim=zlim.x, asp=aspect, axes=axes, xlab=xlab,
                      ylab=ylab, ...)
      graphics::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y, add=TRUE)
    }
  }
  par(oldpar)
  invisible()
}
#' @export
#' @rdname overlay-methods
setGeneric("overlay", function(x, y, ...) standardGeneric("overlay"))
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="nifti", y="nifti"), overlay.nifti)
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="anlz", y="anlz"), overlay.nifti)
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="anlz", y="nifti"), overlay.nifti)
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="nifti", y="anlz"), overlay.nifti)
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="array", y="array"),
          function(x, y, ...) {
            x <- as(x, "nifti")
            y <- as(y, "nifti")
            overlay.nifti(x, y, ...)
          })
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="array", y="nifti"),
          function(x, y, ...) {
            x <- as(x, "nifti")
            overlay.nifti(x, y, ...)
          })
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="nifti", y="array"),
          function(x, y, ...) {
            y <- as(y, "nifti")
            overlay.nifti(x, y, ...)
          })
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="array", y="anlz"),
          function(x, y, ...) {
            x <- as(x, "nifti")
            overlay.nifti(x, y, ...)
          })
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="anlz", y="array"),
          function(x, y, ...) {
            y <- as(y, "nifti")
            overlay.nifti(x, y, ...)
          })
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="afni", y="afni"),
          function(x, y, ...) {
            x <- as(x, "nifti")
            y <- as(y, "nifti")
            overlay.nifti(x, y, ...)
          })
#' @export
#' @rdname overlay-methods
setMethod("overlay", signature(x="afni", y="array"),
          function(x, y, ...) {
            x <- as(x, "nifti")
            overlay.nifti(x, y, ...)
          })

#############################################################################
## orthographic() for class="nifti"
#############################################################################
#' Methods for Function `orthographic' in Package `dcemriS4'
#' 
#' Produce orthographic display for \code{nifti}, \code{anlz} and \code{afni}
#' objects.
#' 
#' 
#' @name orthographic-methods
#' @aliases orthographic orthographic-methods orthographic,nifti-method
#' orthographic,anlz-method orthographic,array-method orthographic,afni-method
#' orthographic.nifti
#' @docType methods
#' @param x is an object of class \code{nifti} or similar.
#' @param y is an object of class \code{nifti} or similar for the overlay.
#' @param xyz is the coordinate for the center of the crosshairs.
#' @param w is the time point to be displayed (4D arrays only).
#' @param col is grayscale (by default).
#' @param col.y is hotmetal (by default).
#' @param zlim is the minimum and maximum \sQuote{z} values passed into
#' \code{image}.
#' @param zlim.y is the minimum and maximum \sQuote{z} values passed into
#' \code{image} for the overlay.
#' @param crosshairs is a logical value for the presence of crosshairs in all
#' three orthogonal planes (default = TRUE).
#' @param col.crosshairs is the color of the crosshairs (default = red).
#' @param xlab is set to "" since all margins are set to zero.
#' @param ylab is set to "" since all margins are set to zero.
#' @param axes is set to \code{FALSE} since all margins are set to zero.
#' @param oma is the size of the outer margins in the \code{par} function.
#' @param mar is the number of lines of margin in the \code{par} function.
#' @param bg is the background color in the \code{par} function.
#' @param text allows the user to specify text to appear in the fourth (unused)
#' pane.
#' @param text.color is the color of the user-specified text (default =
#' \dQuote{white}.
#' @param text.cex is the size of the user-specified text (default = 2).
#' @param \dots other arguments to the \code{image} function may be provided
#' here.
#' @section Methods: \describe{ \item{x = "afni"}{Produce orthographic display
#' for \code{x}.} \item{x = "anlz"}{Produce orthographic display for \code{x}.}
#' \item{x = "array"}{Produce orthographic display for \code{x}.} \item{x =
#' "nifti"}{Produce orthographic display for \code{x}.} }
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{image-methods}}, \code{\link{overlay-methods}}
#' @keywords methods
#' @export
#' @rdname orthographic-methods
orthographic.nifti <- function(x, y=NULL, xyz=NULL, w=1, col=gray(0:64/64),
                               col.y=hotmetal(), zlim=NULL, zlim.y=NULL,
                               crosshairs=TRUE, col.crosshairs="red", 
                               xlab="", ylab="", axes=FALSE,
                               oma=rep(0,4), mar=rep(0,4), bg="black",
                               text=NULL, text.color="white",
                               text.cex=2, ...) {
  if (! is.null(y)) {
    ## both volumes must have the same dimension
    if (! all(dim(x)[1:3] == dim(y)[1:3])) {
      stop("dimensions of \"x\" and \"y\" must be equal")
    }
  }
  X <- nrow(x)
  Y <- ncol(x)
  Z <- nsli(x)
  W <- ntim(x)
  ## Center crosshairs if not specified
  if (is.null(xyz)) {
    xyz <- ceiling(c(X,Y,Z)/2)
  }
  ## check dimensions
  if (X == 0 || Y == 0 || Z == 0) {
    stop("size of NIfTI volume is zero, nothing to plot")
  }
  ## check for z-limits in x; use internal by default
  if (is.null(zlim)) {
    zlim <- c(x@"cal_min", x@"cal_max")
    if (any(!is.finite(zlim)) || diff(zlim) == 0) {
      zlim <- c(x@"glmin", x@"glmax")
    }
    if (any(!is.finite(zlim)) || diff(zlim) == 0) {
      zlim <- range(x, na.rm=TRUE)
    }
  }
  breaks <- c(min(x, zlim, na.rm=TRUE),
              seq(min(zlim, na.rm=TRUE), max(zlim, na.rm=TRUE),
                  length=length(col)-1),
              max(x, zlim, na.rm=TRUE))
  if (! is.null(y) && is.null(zlim.y)) {
    zlim.y <- c(y@"cal_min", y@"cal_max")
    if (max(zlim.y) == 0) {
      zlim.y <- c(x@"glmin", x@"glmax")
    }
  }
  oldpar <- par(no.readonly=TRUE)
  par(mfrow=c(2,2), oma=oma, mar=mar, bg=bg)
  if (is.na(W)) {
    ## Three-dimensional array
    graphics::image(1:X, 1:Z, x[,xyz[2],], col=col, zlim=zlim, breaks=breaks,
                    asp=x@pixdim[4]/x@pixdim[2],
                    xlab=ylab, ylab=xlab, axes=axes, ...)
    if (! is.null(y)) {
      graphics::image(1:X, 1:Z, y[,xyz[2],], col=col.y, zlim=zlim.y, add=TRUE)
    }
    if (crosshairs) {
      abline(h=xyz[3], v=xyz[1], col=col.crosshairs)
    }
    graphics::image(1:Y, 1:Z, x[xyz[1],,], col=col, breaks=breaks,
                    asp=x@pixdim[4]/x@pixdim[3],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (! is.null(y)) {
      graphics::image(1:Y, 1:Z, y[xyz[1],,], col=col.y, zlim=zlim.y, add=TRUE)
    }
    if (crosshairs) {
      abline(h=xyz[3], v=xyz[2], col=col.crosshairs)
    }
    graphics::image(1:X, 1:Y, x[,,xyz[3]], col=col, breaks=breaks,
                    asp=x@pixdim[3]/x@pixdim[2],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (! is.null(y)) {
      graphics::image(1:X, 1:Y, y[,,xyz[3]], col=col.y, zlim=zlim.y, add=TRUE)
    }
    if (crosshairs) {
      abline(h=xyz[2], v=xyz[1], col=col.crosshairs)
    }
  } else {
    ## Four-dimensional array    
    if (w < 1 || w > W) {
      stop("volume \"w\" out of range")
    }
    graphics::image(1:X, 1:Z, x[,xyz[2],,w], col=col, breaks=breaks,
                    asp=x@pixdim[4]/x@pixdim[2],
                    xlab=ylab, ylab=xlab, axes=axes, ...)
    if (! is.null(y)) {
      graphics::image(1:X, 1:Z, y[,xyz[2],], col=col.y, zlim=zlim.y, add=TRUE)
    }
    if (crosshairs) {
      abline(h=xyz[3], v=xyz[1], col=col.crosshairs)
    }
    graphics::image(1:Y, 1:Z, x[xyz[1],,,w], col=col, breaks=breaks,
                    asp=x@pixdim[4]/x@pixdim[3],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (! is.null(y)) {
      graphics::image(1:Y, 1:Z, y[xyz[1],,], col=col.y, zlim=zlim.y, add=TRUE)
    }
    if (crosshairs) {
      abline(h=xyz[3], v=xyz[2], col=col.crosshairs)
    }
    graphics::image(1:X, 1:Y, x[,,xyz[3],w], col=col, breaks=breaks,
                    asp=x@pixdim[3]/x@pixdim[2],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (! is.null(y)) {
      graphics::image(1:X, 1:Y, y[,,xyz[3]], col=col.y, zlim=zlim.y, add=TRUE)
    }
    if (crosshairs) {
      abline(h=xyz[2], v=xyz[1], col=col.crosshairs)
    }
  }
  if (! is.null(text)) {
    ## Add user-supplied text to the "fourth" plot
    graphics::image(1:64, 1:64, matrix(NA, 64, 64), xlab="", ylab="",
                    axes=FALSE)
    text(32, 32, text, col=text.color, cex=text.cex)
  }
  par(oldpar)
  invisible()
}

#' @export
#' @rdname orthographic-methods
setGeneric("orthographic", function(x, ...) standardGeneric("orthographic"))
#' @export
#' @rdname orthographic-methods
setMethod("orthographic", signature(x="nifti"), orthographic.nifti)
#' @export
#' @rdname orthographic-methods
setMethod("orthographic", signature(x="anlz"), orthographic.nifti)
#' @export
#' @rdname orthographic-methods
setMethod("orthographic", signature(x="array"),
          function(x, ...) {
            x <- as(x, "nifti")
            orthographic.nifti(x, ...)
          })
#' @export
#' @rdname orthographic-methods
setMethod("orthographic", signature(x="afni"),
          function(x, ...) {
            x <- as(x, "nifti")
            orthographic.nifti(x, ...)
          })
