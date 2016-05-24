##
## Copyright (c) 2010-2015, Brandon Whitcher
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

#' Read All DICOM Files in a Directory
#' 
#' All DICOM files are imported and a text file summarizing their content
#' recorded.
#' 
#' A \code{for} loop is used to process each DICOM file contained in the
#' directory(ies).  If only a single file is specified in the path,
#' \code{readDICOM} will read that file only.
#' 
#' @aliases readDICOM
#' @param path Path name to the DICOM directory.
#' @param recursive Search recursively down from the given path name.
#' @param exclude Exclude file names containing this character string.
#' @param verbose Flag to provide text-based progress bar.
#' @param counter Ignored.
#' @param ... Arguments to be passed to \code{readDICOMFile}.
#' @return A list structure with two major components: \item{img}{All images
#' associated with the DICOM directory(ies).} \item{hdr}{All header files
#' associated with the DICOM directory(ies).}
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{readDICOMFile}}
#' @references Whitcher, B., V. J. Schmid and A. Thornton (2011).  Working with
#' the DICOM and NIfTI Data Standards in R, \emph{Journal of Statistical
#' Software}, \bold{44} (6), 1--28.  \url{http://www.jstatsoft.org/v44/i06}
#' 
#' Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}
#' @keywords file
#' @examples
#' 
#' ## pixelData = TRUE
#' ## The DICOM image data are read from readDICOM()
#'
#' \dontrun{
#' dcmSphere <- readDICOM(system.file("sphere3", package="oro.dicom"), verbose=TRUE)
#' }
#' 
#' @export readDICOM
readDICOM <- function(path, recursive=TRUE, exclude=NULL, verbose=FALSE,
                      counter, ...) {
  if (length(list.files(path)) == 0 && file.exists(path)) {
    filenames <- path
  } else {
    if (recursive) {
      filenames <- list.files(path, full.names=TRUE, recursive=TRUE)
    } else {
      filenames <- list.files(path, full.names=TRUE)
    }
  }
  if (! is.null(exclude)) {
    filenames <- grep(exclude, filenames, ignore.case=TRUE, value=TRUE,
                      invert=TRUE)
  }
  nfiles <- length(filenames)
  nch <- nchar(as.character(nfiles))
  headers <- images <- vector("list", nfiles)
  names(images) <- names(headers) <- filenames
  if (verbose) {
    cat(" ", nfiles, "files to be processed by readDICOM()", fill=TRUE)
    tpb <- txtProgressBar(min=0, max=nfiles, style=3)
  }
  for (i in 1:nfiles) {
    if (verbose) {
      setTxtProgressBar(tpb, i)
    }
    dcm <- readDICOMFile(filenames[i], ...)
    images[[i]] <- dcm$img
    headers[[i]] <- dcm$hdr
  }
  if (verbose) {
    close(tpb)
  }
  list(hdr=headers, img=images)
}
