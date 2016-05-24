##
## Copyright (c) 2010, Brandon Whitcher
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
##     * Neither the name of Rigorous Analytics Ltd. nor the names of
##       its contributors may be used to endorse or promote products
##       derived from this software without specific prior written
##       permission.
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
## $Id: $
##

#' Convert Direction Cosines to Anatomical Direction
#'
#' For cross-sectional DICOM images the orientation must be derived from the
#' Image Orientation (Patient) direction cosines.
#'
#' C.7.6.2.1.1 Image Position And Image Orientation.  The Image Position
#' (0020,0032) specifies the x, y, and z coordinates of the upper left hand
#' corner of the image; it is the center of the first voxel transmitted.  Image
#' Orientation (0020,0037) specifies the direction cosines of the first row and
#' the first column with respect to the patient.  These Attributes shall be
#' provide as a pair.  Row value for the x, y, and z axes respectively followed
#' by the Column value for the x, y, and z axes respectively.  The direction of
#' the axes is defined fully by the patient's orientation.  The x-axis is
#' increasing to the left hand side of the patient.  The y-axis is increasing
#' to the posterior side of the patient.  The z-axis is increasing toward the
#' head of the patient.  The patient based coordinate system is a right handed
#' system; i.e., the vector cross product of a unit vector along the positive
#' x-axis and a unit vector along the positive y-axis is equal to a unit vector
#' along the positive z-axis.
#'
#' @usage getOrientation(xyz, delta = 0.0001)
#' @param xyz is a vector of direction cosines from
#' \dQuote{ImageOrientationPatient} (0020,0037).
#' @param delta is the tolerance around zero for comparisons.
#' @return Anatomical direction shall be designated by the capital letters:
#' \item{A}{anterior} \item{P}{posterior} \item{R}{right} \item{L}{left}
#' \item{H}{head} \item{F}{foot}
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{swapDimension}}
#' @references \url{http://www.dclunie.com/medical-image-faq/html/part2.html}
#' @export getOrientation
getOrientation <- function(xyz, delta=0.0001) {
  oX <- ifelse(xyz[1] < 0, "R", "L")
  oY <- ifelse(xyz[2] < 0, "A", "P")
  oZ <- ifelse(xyz[3] < 0, "F", "H")
  aX <- abs(xyz[1])
  aY <- abs(xyz[2])
  aZ <- abs(xyz[3])

  orientation <- NULL
  for (i in 1:3) {
    if (aX > delta && aX > aY && aX > aZ) {
      orientation <- paste(orientation, oX, sep="")
      aX <- 0
    } else {
      if (aY > delta && aY > aX && aY > aZ) {
        orientation <- paste(orientation, oY, sep="")
        aY <- 0
      } else {
        if (aZ > delta && aZ > aX && aZ > aY) {
          orientation <- paste(orientation, oZ, sep="")
          aZ <- 0
        }
      }
    }
  }
  return(orientation)
}

#' Reslice Data Volume Using DICOM Header Fields
#'
#' The input data volume (assumed to be three-dimensional) is re-sliced so that
#' each slice is in the axial plane.  Orientation is preserved so that
#' orthographic viewing is standardized.
#'
#' @param img Multidimensional array (assumed to be three-dimensional only).
#' @param dcm DICOM header/image object associated with the multidimensional
#' array.
#' @param digits Number of significant digits used in testing
#' \code{unique}-ness of values in DICOM header fields.
#' @return Multidimensional array with (potentially) permuted dimensions
#' because of the reslicing operation.  An additional attribute
#' \dQuote{\code{pixdim}} is provided in order to facilitate conversion from
#' DICOM to NIFTI/ANALYZE.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{dicom2nifti}}, \code{\link{getOrientation}}
#' @keywords misc
#' @export swapDimension
swapDimension <- function(img, dcm, digits=2) {
  imagePositionPatient <-
    header2matrix(extractHeader(dcm$hdr, "ImagePositionPatient", FALSE), 3)
  if (nrow(imagePositionPatient) != oro.nifti::nsli(img)) {
    imagePositionPatient <- attributes(img)$ipp
  }
  imageOrientationPatient <-
    header2matrix(extractHeader(dcm$hdr, "ImageOrientationPatient", FALSE), 6)
  if (nrow(imageOrientationPatient) != oro.nifti::nsli(img)) {
    imageOrientationPatient <- attributes(img)$iop
  }
  ## Ensure all rows of imageOrientationPatient are identical!
  pixelSpacing <-
    header2matrix(extractHeader(dcm$hdr, "PixelSpacing", FALSE), 2)
  ## Ensure all rows of pixelSpacing are identical!
  sliceThickness <- extractHeader(dcm$hdr, "SliceThickness")
  pixdim <- c(unique(pixelSpacing), unique(sliceThickness))
  ipp.signif <- signif(imageOrientationPatient, digits) # significant digits
  first.row <- getOrientation(unique(signif(ipp.signif))[1:3])
  first.col <- getOrientation(unique(signif(ipp.signif))[4:6])
  if (nchar(first.row) > 1 || nchar(first.col) > 1) {
    warning("Oblique acquisition in ImageOrientationPatient (hope for the best).")
  }
  X <- nrow(img)
  Y <- ncol(img)
  Z <- dim(img)[3]
  W <- dim(img)[4]
  ld <- as.numeric(length(dim(img)))
  ## AXIAL
  if (is.axial(imageOrientationPatient)) {
    if (unlist(strsplit(first.row, ""))[1] %in% c("A","P")) {
      if (ld == 3) {
        index <- c(2,1,3)
      }
      if (ld == 4) {
        index <- c(2,1,3,4)
      }
      img <- aperm(img, index)
      pixdim <- pixdim[index[1:3]]
    }
    if (unlist(strsplit(first.row, ""))[1] == "R") {
      img <- switch(as.character(ld),
                    "3" = img[X:1, , ],
                    "4" = img[X:1, , , ],
                    stop("Dimension \"DIM\" incorrectly specified."))
    }
    if (unlist(strsplit(first.col, ""))[1] == "A") {
      img <- switch(as.character(ld),
                    "3" = img[, Y:1, ],
                    "4" = img[, Y:1, , ])
    }
    ## The z-axis is increasing toward the HEAD of the patient.
    z.index <- order(imagePositionPatient[, 3])
    if (ld == 3) {
      img <- img[, , z.index]
    }
    if (ld == 4) {
      img <- img[, , Z:1, ]
    }
  }
  ## CORONAL
  if (is.coronal(imageOrientationPatient)) {
    if (unlist(strsplit(first.row, ""))[1] %in% c("H","F")) {
      if (ld == 3) {
        index <- c(2,1,3)
      }
      if (ld == 4) {
        index <- c(2,1,3,4)
      }
      img <- aperm(img, index)
      pixdim <- pixdim[index[1:3]]
    }
    if (unlist(strsplit(first.row, ""))[1] == "R") {
      img <- switch(as.character(ld),
                    "3" = img[X:1, , ],
                    "4" = img[X:1, , , ],
                    stop("Dimension \"DIM\" incorrectly specified."))
    }
    if (unlist(strsplit(first.col, ""))[1] == "H") {
      img <- switch(as.character(ld),
                    "3" = img[, Y:1, ],
                    "4" = img[, Y:1, , ])
    }
    ## The y-axis is increasing to the posterior side of the patient.
    z.index <- order(imagePositionPatient[, 2])
    if (ld == 3) {
      img <- img[, , z.index]
      index <- c(1,3,2)
    }
    if (ld == 4) {
      img <- img[, , Z:1, ]
      index <- c(1,3,2,4)
    }
    img <- aperm(img, index) # re-organize orthogonal views
    pixdim <- pixdim[index[1:3]]
  }
  ## SAGITTAL
  if (is.sagittal(imageOrientationPatient)) {
    if (unlist(strsplit(first.row, ""))[1] %in% c("H","F")) {
      if (ld == 3) {
        index <- c(2,1,3)
      }
      if (ld == 4) {
        index <- c(2,1,3,4)
      }
      img <- aperm(img, index)
      pixdim <- pixdim[index[1:3]]
    }
    if (unlist(strsplit(first.row, ""))[1] == "P") {
      img <- switch(as.character(ld),
                    "3" = img[X:1, , ],
                    "4" = img[X:1, , , ],
                    stop("Dimension \"DIM\" incorrectly specified."))
    }
    if (unlist(strsplit(first.col, ""))[1] == "H") {
      img <- switch(as.character(ld),
                    "3" = img[, Y:1, ],
                    "4" = img[, Y:1, , ])
    }
    ## The x-axis is increasing to the left hand side of the patient.
    z.index <- order(imagePositionPatient[, 1])
    if (ld == 3) {
      img <- img[, , z.index]
      index <- c(3,1,2)
    }
    if (ld == 4) {
      img <- img[, , Z:1, ]
      index <- c(3,1,2,4)
    }
    img <- aperm(img, index) # re-organize orthogonal views
    pixdim <- pixdim[index[1:3]]
  }
  imageOrientationPatient <- imageOrientationPatient[z.index, ]
  if (any(is.na(imageOrientationPatient))) {
    stop("Missing values are present in ImageOrientationPatient.")
  }
  imagePositionPatient <- imagePositionPatient[z.index, ]
  if (any(is.na(imagePositionPatient))) {
    stop("Missing values are present in ImagePositionPatient.")
  }
  attr(img, "ipp") <- imagePositionPatient
  attr(img, "iop") <- imageOrientationPatient
  attr(img, "pixdim") <- pixdim
  return(img)
}

#' Orthogonal Planes
#'
#' Functions to test the orientation for a single slice.
#'
#' @aliases is.axial is.coronal is.sagittal
#' @name orthogonal-planes
#' @param imageOrientationPatient A vector of length six taken from the DICOM
#' header field \dQuote{ImageOrientationPatient}.
#' @param axial Characters that are valid in defining an \sQuote{axial} slice.
#' @param coronal Characters that are valid in defining a \sQuote{coronal}
#' slice.
#' @param sagittal Characters that are valid in defining a \sQuote{sagittal}
#' slice.
#' @return Logical value.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{getOrientation}}
#' @keywords misc
#' @examples
#'
#' x <- readDICOMFile(system.file("dcm/Abdo.dcm", package="oro.dicom"))
#' iop <- header2matrix(extractHeader(x$hdr, "ImageOrientationPatient", FALSE), 6)
#' is.axial(iop)
#' is.coronal(iop)
#' is.sagittal(iop)
#'
#' x <- readDICOMFile(system.file("dcm/Spine1.dcm", package="oro.dicom"))
#' iop <- header2matrix(extractHeader(x$hdr, "ImageOrientationPatient", FALSE), 6)
#' is.axial(iop)
#' is.coronal(iop)
#' is.sagittal(iop)
#'
#' @rdname is.axial
#' @export
is.axial <- function(imageOrientationPatient, axial=c("L","R","A","P")) {
  first.row <- getOrientation(imageOrientationPatient[1, 1:3])
  first.col <- getOrientation(imageOrientationPatient[1, 4:6])
  if (nchar(first.row) > 1 || nchar(first.col) > 1) {
    warning("Oblique acquisition in ImageOrientationPatient.")
  }
  return(unlist(strsplit(first.row, ""))[1] %in% axial &
         unlist(strsplit(first.col, ""))[1] %in% axial)
}
#' @rdname is.axial
#' @export
is.coronal <- function(imageOrientationPatient,
                       coronal=c("L","R","H","F")) {
  first.row <- getOrientation(imageOrientationPatient[1, 1:3])
  first.col <- getOrientation(imageOrientationPatient[1, 4:6])
  if (nchar(first.row) > 1 || nchar(first.col) > 1) {
    warning("Oblique acquisition in ImageOrientationPatient.")
  }
  return(unlist(strsplit(first.row, ""))[1] %in% coronal &
         unlist(strsplit(first.col, ""))[1] %in% coronal)
}
#' @rdname is.axial
#' @export
is.sagittal <- function(imageOrientationPatient,
                        sagittal=c("A","P","H","F")) {
  first.row <- getOrientation(imageOrientationPatient[1, 1:3])
  first.col <- getOrientation(imageOrientationPatient[1, 4:6])
  if (nchar(first.row) > 1 || nchar(first.col) > 1) {
    warning("Oblique acquisition in ImageOrientationPatient.")
  }
  return(unlist(strsplit(first.row, ""))[1] %in% sagittal &
         unlist(strsplit(first.col, ""))[1] %in% sagittal)
}
