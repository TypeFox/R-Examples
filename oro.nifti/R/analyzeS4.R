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
## $Id: analyzeS4.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## setClass("anlz")
#############################################################################
#' Class "anlz"
#' 
#' The ANALYZE class for medical imaging data.
#' 
#' @name anlz-class
#' @aliases anlz-class show,anlz-method
#' @param object An object of class \code{anlz}. 
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("anlz", data, dim, dimnames, ...)} or by calling the \code{anlz}
#' function.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\linkS4class{nifti}}, \code{\linkS4class{niftiExtension}}
#' @references ANALYZE 7.5\cr\url{https://rportal.mayo.edu/bir/ANALYZE75.pdf}
#' @keywords classes
#' @examples
#' 
#' showClass("anlz")
#' @section Slots:
#'   \describe{
#'     \item{\code{.Data}:}{Object of class \code{"array"} contains the
#'                          imaging data}
#'     \item{\code{sizeof_hdr}:}{Object of class \code{"numeric"} contains
#'                               the size of the header (= 348)}
#'     \item{\code{data_type}:}{Object of class \code{"character"}}
#'     \item{\code{db_name}:}{Object of class \code{"character"}}
#'     \item{\code{extents}:}{Object of class \code{"numeric"}}
#'     \item{\code{session_error}:}{Object of class \code{"numeric"}}
#'     \item{\code{regular}:}{Object of class \code{"character"}}
#'     \item{\code{hkey_un0}:}{Object of class \code{"character"}}
#'     \item{\code{dim_}:}{Object of class \code{"vector"} contains the
#'                         dimensions of the imaging data}
#'     \item{\code{vox_units}:}{Object of class \code{"character"}}
#'     \item{\code{cal_units}:}{Object of class \code{"character"}}
#'     \item{\code{unused1}:}{Object of class \code{"numeric"}}
#'     \item{\code{datatype}:}{Object of class \code{"numeric"}}
#'     \item{\code{bitpix}:}{Object of class \code{"numeric"} contains the
#'                           number of bits per voxel (pixel)}
#'     \item{\code{dim_un0}:}{Object of class \code{"numeric"}}
#'     \item{\code{pixdim}:}{Object of class \code{"vector"} contains the
#'                           real-world dimensions of the imaging data}
#'     \item{\code{vox_offset}:}{Object of class \code{"numeric"}}
#'     \item{\code{funused1}:}{Object of class \code{"numeric"}}
#'     \item{\code{funused2}:}{Object of class \code{"numeric"}}
#'     \item{\code{funused3}:}{Object of class \code{"numeric"}}
#'     \item{\code{cal_max}:}{Object of class \code{"numeric"} contains the
#'                            maximum display intensity}
#'     \item{\code{cal_min}:}{Object of class \code{"numeric"} contains the
#'                            minimum display intensity}
#'     \item{\code{compressed}:}{Object of class \code{"numeric"}}
#'     \item{\code{verified}:}{Object of class \code{"numeric"}}
#'     \item{\code{glmax}:}{Object of class \code{"numeric"}}
#'     \item{\code{glmin}:}{Object of class \code{"numeric"}}
#'     \item{\code{descrip}:}{Object of class \code{"character"}}
#'     \item{\code{aux_file}:}{Object of class \code{"character"}}
#'     \item{\code{orient}:}{Object of class \code{"character"}}
#'     \item{\code{origin}:}{Object of class \code{"numeric"}}
#'     \item{\code{generated}:}{Object of class \code{"character"}}
#'     \item{\code{scannum}:}{Object of class \code{"character"}}
#'     \item{\code{patient_id}:}{Object of class \code{"character"}}
#'     \item{\code{exp_date}:}{Object of class \code{"character"}}
#'     \item{\code{exp_time}:}{Object of class \code{"character"}}
#'     \item{\code{hist_un0}:}{Object of class \code{"character"}}
#'     \item{\code{views}:}{Object of class \code{"numeric"}}
#'     \item{\code{vols_added}:}{Object of class \code{"numeric"}}
#'     \item{\code{start_field}:}{Object of class \code{"numeric"}}
#'     \item{\code{field_skip}:}{Object of class \code{"numeric"}}
#'     \item{\code{omax}:}{Object of class \code{"numeric"}}
#'     \item{\code{omin}:}{Object of class \code{"numeric"}}
#'     \item{\code{smax}:}{Object of class \code{"numeric"}}
#'     \item{\code{smin}:}{Object of class \code{"numeric"}}
#'   }
#' @section Extends:
#'   Class \code{"\linkS4class{array}"}, from data part.\cr
#'   Class \code{"\linkS4class{matrix}"}, by class \dQuote{array}, distance 2,
#'   with explicit test and coerce. \cr
#'   Class \code{"\linkS4class{structure}"}, by class \dQuote{array}, distance 2.\cr
#'   Class \code{"\linkS4class{vector}"}, by class \dQuote{array}, distance 3,
#'   with explicit coerce.\cr
#'   Class \code{"\linkS4class{vector}"}, by class \dQuote{array}, distance 5,
#'   with explicit test and coerce.
#' @section Methods:
#'   \describe{
#'     \item{image}{\code{signature(x = "anlz")}: diplays the image(s).}
#'     \item{show}{\code{signature(object = "anlz")}: prints out a summary
#'                 of the imaging data.}
#'   }
#' @rdname anlz-class
#' @export 
setClass("anlz",
         representation("sizeof_hdr" = "numeric",
                        "data_type" = "character",
                        "db_name" = "character",
                        "extents" = "numeric",
                        "session_error" = "numeric",
                        "regular" = "character",
                        "hkey_un0" = "character",
                        "dim_" = "vector",
                        "vox_units" = "character",
                        "cal_units" = "character",
                        "unused1" = "numeric",
                        "datatype" = "numeric",
                        "bitpix" = "numeric",
                        "dim_un0" = "numeric",
                        "pixdim" = "vector",
                        "vox_offset" = "numeric",
                        "funused1" = "numeric",
                        "funused2" = "numeric",
                        "funused3" = "numeric",
                        "cal_max" = "numeric",
                        "cal_min" = "numeric",
                        "compressed" = "numeric",
                        "verified" = "numeric",
                        "glmax" = "numeric",
                        "glmin" = "numeric",
                        "descrip" = "character",
                        "aux_file" = "character",
                        "orient" = "character",
                        "origin" = "numeric",
                        "generated" = "character",
                        "scannum" = "character",
                        "patient_id" = "character",
                        "exp_date" = "character",
                        "exp_time" = "character",
                        "hist_un0" = "character",
                        "views" = "numeric",
                        "vols_added" = "numeric",
                        "start_field" = "numeric",
                        "field_skip" = "numeric",
                        "omax" = "numeric",
                        "omin" = "numeric",
                        "smax" = "numeric",
                        "smin" = "numeric"),
         prototype("sizeof_hdr" = 348,
                   "data_type" = "",
                   "db_name" = "",
                   "extents" = numeric(1),
                   "session_error" = numeric(1),
                   "regular" = "r",
                   "hkey_un0" = "",
                   "dim_" = numeric(8),
                   "vox_units" = "mm",
                   "cal_units" = "",
                   "unused1" = numeric(1),
                   "datatype" = 2,
                   "bitpix" = 8,
                   "dim_un0" = numeric(1),
                   "pixdim" = numeric(8),
                   "vox_offset" = numeric(1),
                   "funused1" = numeric(1),
                   "funused2" = numeric(1),
                   "funused3" = numeric(1),
                   "cal_max" = numeric(1),
                   "cal_min" = numeric(1),
                   "compressed" = numeric(1),
                   "verified" = numeric(1),
                   "glmax" = numeric(1),
                   "glmin" = numeric(1),
                   "descrip" = "",
                   "aux_file" = "",
                   "orient" = "0",
                   "origin" = numeric(5),
                   "generated" = "",
                   "scannum" = "",
                   "patient_id" = "",
                   "exp_date" = "",
                   "exp_time" = "",
                   "hist_un0" = "",
                   "views" = numeric(1),
                   "vols_added" = numeric(1),
                   "start_field" = numeric(1),
                   "field_skip" = numeric(1),
                   "omax" = numeric(1),
                   "omin" = numeric(1),
                   "smax" = numeric(1),
                   "smin" = numeric(1)),
         contains="array")

#############################################################################
## setMethod("show", "anlz")
#############################################################################
#' @aliases show,anlz-method
#' @rdname anlz-class
setMethod("show", "anlz", function(object) {
  cat("ANALYZE 7.5 format", fill=TRUE)
  cat("  Type            :", class(object), fill=TRUE)
  cat("  Data Type       : ", object@"datatype", " (",
      convert.datatype.anlz(object@"datatype"), ")", sep="", fill=TRUE)
  cat("  Bits per Pixel  :", object@"bitpix", fill=TRUE)
  cat("  Orient          : ", object@"orient", " (",
      convert.orient.anlz(object@"orient"), ")", sep="", fill=TRUE)
  cat("  Dimension       :",
      paste(object@"dim_"[2:(1+object@"dim_"[1])], collapse=" x "), fill=TRUE)
  cat("  Pixel Dimension :",
      paste(round(object@"pixdim"[2:(1+object@"dim_"[1])], 2),
            collapse=" x "), fill=TRUE)
  cat("  Voxel Units     :", object@"vox_units", fill=TRUE)
})

#############################################################################
## setValidity("anlz")
#############################################################################

setValidity("anlz", function(object) {
  retval <- NULL
  indices <- 2:(1+object@"dim_"[1])
  ## sizeof_hdr must be 348
  if (object@"sizeof_hdr" != 348) {
    retval <- c(retval, "sizeof_hdr != 348")
  }
  ## datatype needed to specify type of image data
  if (!object@"datatype" %in% c(0, 2^(0:7), 255)) {
    retval <- c(retval, "datatype not recognized")
  }
  ## bitpix should correspond correctly to datatype
  ## 
  ## dim should be non-zero for dim[1] dimensions
  if (!all(as.logical(object@"dim_"[indices]))) {
    retval <- c(retval, "dim[1]/dim mismatch")
  }
  ## number of data dimensions should match dim[1]
  if (length(indices) != length(dim(object@.Data))) {
    retval <- c(retval, "dim[1]/img mismatch")
  }
  ## pixdim[n] required when dim[n] is required
  if (!all(as.logical(object@"dim_"[indices]) &&
           as.logical(object@"pixdim"[indices]))) {
    retval <- c(retval, "dim/pixdim mismatch")
  }
  ## data dimensions should match dim 
  if (!all(object@"dim_"[indices] == dim(object@.Data))) {
    retval <- c(retval, "dim/img mismatch")
  }
  if (is.null(retval)) {
    return(TRUE)
  } else {
    return(retval)
  }
})

#############################################################################
## anlz()
#############################################################################
#' @title Constructor for Analyze
#' 
#' @description Constructor for Analyze class objects.
#' 
#' @aliases anlz
#' @param img is a multidimensional array of data.
#' @param dim is the dimension of the data (default = \code{missing}).
#' @param datatype is an integer that denotes the type of data contained in
#' each voxel.  See the function \code{convert.datatype.anlz} or the 
#' ANALYZE documentation for more details.
#' @param \dots allows for additional \sQuote{slots} to be specified.
#' @return An object of class \code{anlz}.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\linkS4class{anlz}}, \code{\linkS4class{nifti}},
#' \code{\link{nifti}}, \code{\link{convert.datatype.anlz}}
#' @references ANALYZE 7.5\cr \url{https://rportal.mayo.edu/bir/ANALYZE75.pdf}
#' @examples
#' 
#' aim <- anlz() # default
#' 
#' @export anlz
anlz <- function(img=array(0, dim=rep(1,4)), dim, datatype=2, ...) {
  if (missing(dim)) {
    if (is.array(img)) {
      dim <- base::dim(img)
    } else {
      dim <- c(1, length(img))
    }
  }
  ld <- length(dim)
  if (ld < 3) {
    stop(sprintf("length(dim) must be at least 3 and is %d.", ld))
  }
  x <- c(length(dim), dim[1], dim[2], dim[3],
         ifelse(is.na(dim[4]), 1, dim[4]), rep(1,3))
  y <- c(0.0, rep(1.0, length(dim)), rep(0.0, 7 - length(dim)))
  cal.max <- max(img, na.rm=TRUE)
  cal.min <- min(img, na.rm=TRUE)
  ## Set datatype
  switch(as.character(datatype),
         "1" = bitpix <- 1,
         "2" = bitpix <- 8,
         "4" = bitpix <- 16,
         "8" = bitpix <- 32,
         "16" = bitpix <- 32,
         "32" = bitpix <- 64,
         "64" = bitpix <- 64,
         "512" = bitpix <- 16,
         stop(paste("Data type", datatype, "unsupported."))
         )
  obj <- new("anlz", .Data=array(img, dim=dim), "dim_"=x, "pixdim"=y,
             "cal_max"=cal.max, "cal_min"=cal.min, "datatype"=datatype,
             "bitpix"=bitpix, ...)
  validANALYZE <- getValidity(getClassDef("anlz"))
  validANALYZE(obj)
  return(obj)
}
