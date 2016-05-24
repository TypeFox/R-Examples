## Copyright (c) 2016, James P. Howard, II <jh@jameshoward.us>
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
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

#' @name resizeImage
#' @rdname resizeImage
#'
#' @title Image resizing
#'
#' @description
#' Resize images using nearest neighbor and
#'
#' @param imx a 3-dimensional array containing image data
#' @param width the new width
#' @param height the new height
#'
#' @details
#' The \var{resizeImageNN} function uses the nearest neighbor method to
#' resize the image.  Also, \var{resizeImageBL} uses bilinear
#' interpolation to resize the image.
#'
#' @return a three-dimensional array containing the resized image.
#'
#' @family interpolation

#' @export
resizeImageNN <- function(imx, width, height) {
    imx.dim <- dim(imx)
    layers <- imx.dim[3]
    w.scale <- width / (imx.dim[1] - 1)
    h.scale <- height / (imx.dim[2] - 1)

    ## Create a new array to store the new image
    nn <- array(0, c(height, width, layers))

    for(l in 1:layers)
        for(h in 0:height) {
            y <- round(h / h.scale) + 1
            for(w in 0:width) {
                x <- round(w / w.scale) + 1
                nn[h, w, l] <- imx[y, x, l]
            }
        }

    return(nn)
}

#' @rdname resizeImage
#' @export
resizeImageBL <- function(imx, width, height) {
    imx.dim <- dim(imx)
    layers <- imx.dim[3]
    w.orig <- imx.dim[1]
    h.orig <- imx.dim[2]

    w.scale <- (w.orig - 1) / width
    h.scale <- (h.orig - 1) / height

    ## Create a new array to store the new image
    bl <- array(0, c(height, width, layers))

    zb <- matrix(0, 2, 2)
    for(l in 1:layers)
        for(h in 1:height) {
            y <- floor(h.scale * h)
            if(y == 0) y <- 1
            if(y == h.orig) y <- h.orig - 1

            yb <- c(y, y + 1)
            yl <- h.scale * h

            for(w in 1:width) {
                x <- floor(w.scale * w)
                if(x == 0) x <- 1
                if(x == w.orig) x <- w.orig - 1

                xb <- c(x, x + 1)
                xl <- w.scale * w

                zb[1, 1] <- imx[y, x, l]
                zb[1, 2] <- imx[y, x + 1, l]
                zb[2, 1] <- imx[y + 1, x, l]
                zb[2, 2] <- imx[y + 1, x + 1, l]

                bl[h, w, l] <- bilinear(xb, yb, zb, xl, yl)
            }
        }

    return(bl)
}
