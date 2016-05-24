##
## Copyright (c) 2010-2014 Brandon Whitcher
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

#' Create Arrays from DICOM Headers/Images
#' 
#' A DICOM list structure is used to produce a multi-dimensional array
#' representing a single acquisition of medical imaging data.
#' 
#' 
#' @aliases create3D create4D
#' @param dcm is the DICOM list structure (if \code{pixelData} = \code{TRUE})
#' or the DICOM header information (if \code{pixelData} = \code{FALSE}).
#' @param mode is a valid character string for \code{storage.mode}.
#' @param transpose is available in order to switch the definition of rows and
#' columns from DICOM (default = \code{TRUE}.
#' @param pixelData is a logical variable (default = \code{TRUE}) that is
#' associated with the DICOM image data being pre-loaded.
#' @param mosaic is a logical variable (default = \code{FALSE}) to denote
#' storage of the data in Siemens \sQuote{Mosaic} format.
#' @param mosaicXY is a vector of length two that provides the (x,y) dimensions
#' of the individual images.  Default behavior is to use the AcquisitonMatrix
#' to determine the (x,y) values.
#' @param sequence is a logical variable (default = \code{FALSE}) on whether to
#' look in SequenceItem entries for DICOM header information.
#' @param nslices is the third dimension of the array.  Attempts are made to
#' determine this number from the DICOM data.
#' @param ntimes is the fourth dimension of the array.  Attempts are made to
#' determine this number from the DICOM data.
#' @param instance is a logical variable (default = \code{TRUE}) that
#' determines whether or not to access the \code{InstanceNumber} field in the
#' DICOM header to help order the slices.
#' @return Multi-dimensional array of medical imaging data.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{array}}, \code{\link{readDICOM}},
#' \code{\link{storage.mode}}
#' @references Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}
#' @examples
#' 
#' load(system.file("hk-40/hk40.RData", package="oro.dicom"))
#' dcmList <- hk40
#' dcmImage <- create3D(dcmList)
#' image(dcmImage[,,1], col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
#'       main=paste("First Slice from HK-40"))
#' imagePositionPatient <- attributes(dcmImage)$ipp
#' dSL <- abs(diff(imagePositionPatient[,3]))
#' plot(dSL, ylim=range(range(dSL) * 1.5, 0, 10), xlab="Image", ylab="mm",
#'      main="Difference in Slice Location")
#' 
#' \dontrun{
#' ## pixelData = FALSE
#' ## The DICOM image data are read from create3D()
#' ## This may save on memory for large batches of DICOM data
#' dcmList <- readDICOM(system.file("hk-40", package="oro.dicom"),
#'                      pixelData=FALSE)
#' dcmImage <- create3D(dcmList, pixelData=FALSE)
#' image(dcmImage[,,1], col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
#'       main=paste("First Slice from HK-40 (again)"))
#' }
#' ## mosaic = TRUE
#' mosaicFile <- system.file("dcm/MR-sonata-3D-as-Tile.dcm", package="oro.dicom")
#' dcm <- readDICOMFile(mosaicFile)
#' image(t(dcm$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
#'       main="Siemens MOSAIC")
#' dcmImage <- create3D(dcm, mode="integer", mosaic=TRUE)
#' z <- trunc(dim(dcmImage)[3]/2)
#' image(dcmImage[,,z], col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
#'       main=paste("Slice", z, "from Siemens MOSAIC"))
#' 
#' @export create3D
create3D <- function(dcm, mode="integer", transpose=TRUE, pixelData=TRUE,
                     mosaic=FALSE, mosaicXY=NULL, sequence=FALSE) {
  if (pixelData) {
    if (is.null(dcm$hdr)) {
      stop("DICOM \"hdr\" information is not present.")
    }
    if (is.null(dcm$img)) {
      stop("DICOM \"img\" information is not present.")
    }
  } else {
    if (is.null(dcm$img)) {
      dcm <- list(hdr=dcm, img=NULL) # Only a list of headers as input
    }
  }
  X <- unique(extractHeader(dcm$hdr, "Rows", inSequence=sequence))
  if (length(X) != 1) {
    stop("Row lengths are not identical.")
  }
  Y <- unique(extractHeader(dcm$hdr, "Columns", inSequence=sequence))
  if (length(Y) != 1) {
    stop("Column lengths are not identical.")
  }
  if (mosaic) {
    if (is.null(dim(dcm$img))) {
      stop("Multiple MOSAIC files detected, please use create4D()")
    }
    if (is.null(mosaicXY)) {
      acquisitionMatrix <-
        header2matrix(extractHeader(dcm$hdr, "AcquisitionMatrix", FALSE), 4)
      x <- acquisitionMatrix[1,1]
      y <- acquisitionMatrix[1,4]
      if (is.na(x) || is.na(y)) {
        stop("Missing AcquisitionMatrix, please specify \"mosaicXY\"")
      }
    } else {
      x <- mosaicXY[1]
      y <- mosaicXY[2]
    }
    z <- (X / x) * (Y / y)
    img <- array(0, c(x,y,z))
    k <- 1
    for (i in (X/x):1) {
      for (j in 1:(Y/y)) {
        img[,,k] <- dcm$img[((i-1) * x) + 1:x, ((j-1) * y) + 1:y]
        k <- k + 1
      }
    }
    storage.mode(img) <- mode
    imagePositionPatient <- cbind(X, Y, 1:z)
    imageOrientationPatient <-
      header2matrix(extractHeader(dcm$hdr, "ImageOrientationPatient", FALSE), 6)
    imageOrientationPatient <- matrix(imageOrientationPatient, z, 6, byrow=TRUE)
  } else {
    ## Check if the DICOM list has length > 1
    Z <- ifelse(is.null(dim(dcm$img)), length(dcm$hdr), 1)
    img <- array(0, c(X,Y,Z))
    storage.mode(img) <- mode
    imagePositionPatient <-
      header2matrix(extractHeader(dcm$hdr, "ImagePositionPatient", FALSE), 3)
    if (any(is.na(imagePositionPatient))) {
      stop("Missing values detected in ImagePositionPatient.")
    }
    imageOrientationPatient <-
      header2matrix(extractHeader(dcm$hdr, "ImageOrientationPatient", FALSE), 6)
    if (any(is.na(imageOrientationPatient))) {
      stop("Missing values detected in ImageOrientationPatient.")
    }
      movingDimensions <- apply(imagePositionPatient, 2,
                              function(j) any(diff(j) != 0))
    if (sum(movingDimensions) != 1) {
      warning("ImagePositionPatient is moving in more than one dimension.")
    }
    if (pixelData) {
        for (z in 1:Z) {
            img[,,z] <- dcm$img[[z]]
        }
    } else {
        for (z in 1:Z) {
            img[,,z] <- readDICOMFile(names(dcm$hdr)[z])$img
        }
    }
  }
  if (transpose) {
    img <- aperm(img, c(2,1,3))
  }
  attr(img, "ipp") <- imagePositionPatient
  attr(img, "iop") <- imageOrientationPatient
  return(img)
}
#' @rdname create3D
#' @export create4D
create4D <- function(dcm, mode="integer", transpose=TRUE, pixelData=TRUE,
                     mosaic=FALSE, mosaicXY=NULL, nslices=NULL,
                     ntimes=NULL, instance=TRUE, sequence=FALSE) {
  if (pixelData) {
    if (is.null(dcm$hdr)) {
      stop("DICOM \"hdr\" information is not present.")
    }
    if (is.null(dcm$img)) {
      stop("DICOM \"img\" information is not present.")
    }
  } else {
    if (is.null(dcm$img)) {
      dcm <- list(hdr=dcm, img=NULL) # Only a list of headers as input
    }
  }
  X <- unique(extractHeader(dcm$hdr, "Rows", inSequence=sequence))
  if (length(X) != 1) {
    stop("Row lengths are not identical.")
  }
  Y <- unique(extractHeader(dcm$hdr, "Columns", inSequence=sequence))
  if (length(Y) != 1) {
    stop("Column lengths are not identical.")
  }
  ## Check if the DICOM list has length > 1
  imagePositionPatient <-
    header2matrix(extractHeader(dcm$hdr, "ImagePositionPatient", FALSE), 3)
  if (any(is.na(imagePositionPatient))) {
    stop("Missing values detected in ImagePositionPatient.")
  }
  if (mosaic) {
    if (is.null(mosaicXY)) {
      acquisitionMatrix <-
        header2matrix(extractHeader(dcm$hdr, "AcquisitionMatrix", FALSE), 4)
      x <- acquisitionMatrix[1,1]
      y <- acquisitionMatrix[1,4]
      if (is.na(x) || x == 0 || is.na(y) || y == 0) {
        stop("Missing AcquisitionMatrix, please specify \"mosaicXY\".")
      }
      if (! all(trunc(X/x) == X/x, trunc(Y/y) == Y/y)) {
        if (! all(trunc(X/y) == X/y, trunc(Y/x) == Y/x)) {
          y <- acquisitionMatrix[1,1]
          x <- acquisitionMatrix[1,4]
          warning("AcquisitionMatrix has been transposed.")
        } else {
          stop("AcquisitionMatrix does not make sense, please specify \"mosaicXY\".")
        }
      }
    } else {
      x <- mosaicXY[1]
      y <- mosaicXY[2]
    }
    z <- (X/x) * (Y/y)
    ntimes <- length(dcm$hdr)
    img <- array(0, c(x,y,z,ntimes))
    storage.mode(img) <- mode
    if (pixelData) {
      for (w in 1:ntimes) {
        k <- 1
        for (i in (X/x):1) {
          for (j in 1:(Y/y)) {
            img[,,k,w] <- dcm$img[[w]][((i-1)*x)+1:x, ((j-1)*y)+1:y]
            k <- k+1
          }
        }
      }
    } else {
      for (w in 1:ntimes) {
        k <- 1
        dicomInfoImage <- readDICOMFile(names(dcm$hdr)[w])$img
        for (i in (X/x):1) {
          for (j in 1:(Y/y)) {
            img[,,k,w] <- dicomInfoImage[((i-1)*x)+1:x, ((j-1)*y)+1:y]
            k <- k+1
          }
        }
      }
    }
  } else {
    Z <- ifelse(is.null(dim(dcm$img)), length(dcm$hdr), 1)
    if (Z == 1) {
      warning(paste("Number of DICOM images is", Z, ".", sep=""))
    }
    ## cat("## X =", X, "Y =", Y, "Z =", Z, fill=TRUE)
    movingDimensions <- apply(imagePositionPatient, 2,
                              function(j) any(diff(j) != 0))
    if (sum(movingDimensions) > 1) {
      warning("ImagePositionPatient indicates oblique slices, assuming transverse acquisition.")
      movingDimensions <- c(FALSE, FALSE, TRUE)
    }
    ## Guess number of slices
    if (is.null(nslices)) {
      nslices <- length(unique(imagePositionPatient[,movingDimensions]))
    }
    ## cat("## nslices =", nslices, fill=TRUE)
    if (is.null(nslices)) {
      stop("The number of slices has not been specified/determined.")
    }
    ## Guess the slice order
    instanceNumber <- extractHeader(dcm$hdr, "InstanceNumber")
    if (instance && length(unique(instanceNumber)) == length(dcm$hdr)) {
      index <- order(instanceNumber)
    } else {
      warning("No unique slice ordering found in InstanceNumber.")
      index <- 1:Z
    }
    img <- array(0, c(X,Y,nslices,Z/nslices))
    storage.mode(img) <- mode
    if (pixelData) {
      for (z in 1:Z) {
        zz <- (z - 1) %% nslices + 1
        ww <- (z - 1) %/% nslices + 1
        img[,,zz,ww] <- dcm$img[[index[z]]]
      }
    } else {
      for (z in 1:Z) {
        zz <- (z - 1) %% nslices + 1
        ww <- (z - 1) %/% nslices + 1
        img[,,zz,ww] <- readDICOMFile(names(dcm$hdr)[index[z]])$img
      }
    }
  }
  if (transpose) {
    img <- aperm(img, c(2,1,3,4))
  }
  attr(img, "ipp") <- imagePositionPatient
  return(img)
}
