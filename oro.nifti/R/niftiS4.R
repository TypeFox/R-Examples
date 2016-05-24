##
## Copyright (c) 2009-2011 Brandon Whitcher and Volker Schmid
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
## $Id: niftiS4.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## setClass("niftiExtensionSection")
#############################################################################
#' @title Class "niftiExtensionSection"
#'
#' @description A \code{niftiExtensionSection} contains the fields that conform
#'   to the NIfTI standard regarding header extensions.  A \code{niftiExtension}
#'   is composed of one or more of these objects.
#'
#' @name niftiExtensionSection-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("niftiExtensionSection", data, dim, dimnames, ...)}.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com},\cr Andrew Thornton
#' \email{zeripath@@users.sourcefore.net}
#' @seealso \code{\linkS4class{niftiExtension}}, \code{\linkS4class{nifti}}
#' @references NIfTI-1\cr \url{http://nifti.nimh.nih.gov/}
#' @keywords classes
#' @examples
#'
#' showClass("niftiExtensionSection")
#' @export
setClass("niftiExtensionSection",
         representation(esize="numeric",
                        ecode="numeric",
                        edata="character"),
         prototype(esize=numeric(1),
                   ecode=numeric(1),
                   edata=""))

#############################################################################
## setClass("nifti")
#############################################################################
#' @name nifti-class
#' @title Class "nifti"
#'
#' @description The NIfTI class for medical imaging data.
#'
#' @aliases nifti-class show,nifti-method
#' @param object An object of class \code{nifti}.
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("nifti", data, dim, dimnames, ...)} or by calling the \code{nifti}
#' function.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com},\cr Andrew Thornton
#' \email{zeripath@@users.sourcefore.net}
#' @seealso \code{\linkS4class{anlz}}, \code{\linkS4class{niftiExtension}},
#' \code{\linkS4class{niftiAuditTrail}}
#' @references NIfTI-1\cr \url{http://nifti.nimh.nih.gov/}
#' @keywords classes
#' @examples
#'
#' showClass("nifti")
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
#'     \item{\code{dim_info}:}{Object of class \code{"numeric"} contains
#'                             MRI slice ordering}
#'     \item{\code{dim_}:}{Object of class \code{"vector"} contains the
#'                         dimensions of the imaging data}
#'     \item{\code{intent_p1}:}{Object of class \code{"numeric"}}
#'     \item{\code{intent_p2}:}{Object of class \code{"numeric"}}
#'     \item{\code{intent_p3}:}{Object of class \code{"numeric"}}
#'     \item{\code{intent_code}:}{Object of class \code{"numeric"}}
#'     \item{\code{datatype}:}{Object of class \code{"numeric"}}
#'     \item{\code{bitpix}:}{Object of class \code{"numeric"} contains the
#'                           number of bits per voxel (pixel)}
#'     \item{\code{slice_start}:}{Object of class \code{"numeric"}}
#'     \item{\code{pixdim}:}{Object of class \code{"vector"} contains the
#'                           real-world dimensions of the imaging data}
#'     \item{\code{vox_offset}:}{Object of class \code{"numeric"} contains
#'                               the voxel offset (= 352 when no extensions exist)}
#'     \item{\code{scl_slope}:}{Object of class \code{"numeric"}}
#'     \item{\code{scl_inter}:}{Object of class \code{"numeric"}}
#'     \item{\code{slice_end}:}{Object of class \code{"numeric"}}
#'     \item{\code{slice_code}:}{Object of class \code{"numeric"}}
#'     \item{\code{xyzt_units}:}{Object of class \code{"numeric"}}
#'     \item{\code{cal_max}:}{Object of class \code{"numeric"} contains the
#'                            maximum display intensity}
#'     \item{\code{cal_min}:}{Object of class \code{"numeric"} contains the
#'                            minimum display intensity}
#'     \item{\code{slice_duration}:}{Object of class \code{"numeric"}}
#'     \item{\code{toffset}:}{Object of class \code{"numeric"}}
#'     \item{\code{glmax}:}{Object of class \code{"numeric"}}
#'     \item{\code{glmin}:}{Object of class \code{"numeric"}}
#'     \item{\code{descrip}:}{Object of class \code{"character"}}
#'     \item{\code{aux_file}:}{Object of class \code{"character"}}
#'     \item{\code{qform_code}:}{Object of class \code{"numeric"}}
#'     \item{\code{sform_code}:}{Object of class \code{"numeric"}}
#'     \item{\code{quatern_b}:}{Object of class \code{"numeric"}}
#'     \item{\code{quatern_c}:}{Object of class \code{"numeric"}}
#'     \item{\code{quatern_d}:}{Object of class \code{"numeric"}}
#'     \item{\code{qoffset_x}:}{Object of class \code{"numeric"}}
#'     \item{\code{qoffset_y}:}{Object of class \code{"numeric"}}
#'     \item{\code{qoffset_z}:}{Object of class \code{"numeric"}}
#'     \item{\code{srow_x}:}{Object of class \code{"vector"}}
#'     \item{\code{srow_y}:}{Object of class \code{"vector"}}
#'     \item{\code{srow_z}:}{Object of class \code{"vector"}}
#'     \item{\code{intent_name}:}{Object of class \code{"character"}}
#'     \item{\code{magic}:}{Object of class \code{"character"}}
#'     \item{\code{extender}:}{Object of class \code{"vector"}}
#'     \item{\code{reoriented}:}{Object of class \code{"logical"}}
#'   }
#' @section Extends:
#'   Class \code{"\linkS4class{array}"}, from data part.\cr
#'   Class \code{"\linkS4class{matrix}"}, by class \dQuote{array}, distance 2,
#'   with explicit test and coerce.\cr
#'   Class \code{"\linkS4class{structure}"}, by class \dQuote{array}, distance
#'   2.\cr
#'   Class \code{"\linkS4class{vector}"}, by class \dQuote{array}, distance 3,
#'   with explicit coerce.\cr
#'   Class \code{"\linkS4class{vector}"}, by class \dQuote{array}, distance 5,
#'   with explicit test and coerce.
#'
#' @section Methods:
#'   \describe{
#'     \item{image}{\code{signature(x = "nifti")}: diplays the image(s).}
#'     \item{orthographic}{\code{signature(x = "nifti")}: displays the image(s).}
#'     \item{overlay}{\code{signature(x = "nifti", y = "nifti")}: displays
#'                    the image(s).}
#'     \item{show}{\code{signature(object = "nifti")}: prints out a summary
#'                 of the imaging data.}
#'   }
#' @export
setClass("nifti",
         representation("sizeof_hdr"="numeric",
                        "data_type"="character",
                        "db_name"="character",
                        "extents"="numeric",
                        "session_error"="numeric",
                        "regular"="character",
                        "dim_info"="character",
                        "dim_"="vector",
                        "intent_p1"="numeric",
                        "intent_p2"="numeric",
                        "intent_p3"="numeric",
                        "intent_code"="numeric",
                        "datatype"="numeric",
                        "bitpix"="numeric",
                        "slice_start"="numeric",
                        "pixdim"="vector",
                        "vox_offset"="numeric",
                        "scl_slope"="numeric",
                        "scl_inter"="numeric",
                        "slice_end"="numeric",
                        "slice_code"="numeric", # character?
                        "xyzt_units"="numeric", # character?
                        "cal_max"="numeric",
                        "cal_min"="numeric",
                        "slice_duration"="numeric",
                        "toffset"="numeric",
                        "glmax"="numeric",
                        "glmin"="numeric",
                        "descrip"="character",
                        "aux_file"="character",
                        "qform_code"="numeric",
                        "sform_code"="numeric",
                        "quatern_b"="numeric",
                        "quatern_c"="numeric",
                        "quatern_d"="numeric",
                        "qoffset_x"="numeric",
                        "qoffset_y"="numeric",
                        "qoffset_z"="numeric",
                        "srow_x"="vector",
                        "srow_y"="vector",
                        "srow_z"="vector",
                        "intent_name"="character",
                        "magic"="character",
                        "extender"="vector",
                        "reoriented"="logical"),
         prototype("sizeof_hdr"=348,
                   "data_type"="",
                   "db_name"="",
                   "extents"=numeric(1),
                   "session_error"=numeric(1),
                   "regular"="",
                   "dim_info"="",
                   "dim_"=numeric(8),
                   "intent_p1"=numeric(1),
                   "intent_p2"=numeric(1),
                   "intent_p3"=numeric(1),
                   "intent_code"=numeric(1),
                   "datatype"=2,
                   "bitpix"=8,
                   "slice_start"=numeric(1),
                   "pixdim"=numeric(8),
                   "vox_offset"=352,
                   "scl_slope"=numeric(1),
                   "scl_inter"=numeric(1),
                   "slice_end"=numeric(1),
                   "slice_code"=numeric(1),
                   "xyzt_units"=numeric(1),
                   "cal_max"=numeric(1),
                   "cal_min"=numeric(1),
                   "slice_duration"=numeric(1),
                   "toffset"=numeric(1),
                   "glmax"=numeric(1),
                   "glmin"=numeric(1),
                   "descrip"="",
                   "aux_file"="",
                   "qform_code"=numeric(1),
                   "sform_code"=numeric(1),
                   "quatern_b"=numeric(1),
                   "quatern_c"=numeric(1),
                   "quatern_d"=numeric(1),
                   "qoffset_x"=numeric(1),
                   "qoffset_y"=numeric(1),
                   "qoffset_z"=numeric(1),
                   "srow_x"=numeric(4),
                   "srow_y"=numeric(4),
                   "srow_z"=numeric(4),
                   "intent_name"="",
                   "magic"="n+1",
                   "extender"=numeric(4),
                   "reoriented"=FALSE),
         contains="array")

#############################################################################
## setMethod("show", "nifti")
#############################################################################
#' @aliases show,nifti-method
#' @rdname nifti-class
setMethod("show", "nifti", function(object) {
  cat("NIfTI-1 format", fill=TRUE)
  cat("  Type            :", class(object), fill=TRUE)
  cat("  Data Type       : ", object@"datatype",
      " (", convert.datatype(object@datatype), ")", sep="", fill=TRUE)
  cat("  Bits per Pixel  :", object@bitpix, fill=TRUE)
  cat("  Slice Code      : ", object@"slice_code",
      " (", convert.slice(object@"slice_code"), ")", sep="", fill=TRUE)
  cat("  Intent Code     : ", object@"intent_code",
      " (", convert.intent(object@"intent_code"), ")", sep="", fill=TRUE)
  cat("  Qform Code      : ", object@"qform_code",
      " (", convert.form(object@"qform_code"), ")", sep="", fill=TRUE)
  cat("  Sform Code      : ", object@"sform_code",
      " (", convert.form(object@"sform_code"), ")", sep="", fill=TRUE)
  cat("  Dimension       :",
      paste(object@"dim_"[2:(1+object@"dim_"[1])], collapse=" x "),
      fill=TRUE)
  cat("  Pixel Dimension :",
      paste(round(pixdim(object)[2:(1+object@"dim_"[1])],2), collapse=" x "),
      fill=TRUE)
  cat("  Voxel Units     :", convert.units(xyzt2space(object@"xyzt_units")),
      fill=TRUE)
  cat("  Time Units      :", convert.units(xyzt2time(object@"xyzt_units")),
      fill=TRUE)
})

#############################################################################
## setValidity("nifti")
#############################################################################
setValidity("nifti", function(object) {
  retval <- NULL
  indices <- 1 + 1:object@"dim_"[1]
  ## sizeof_hdr must be 348
  if (object@"sizeof_hdr" != 348) {
    retval <- c(retval, "sizeof_hdr != 348\n")
  }
  ## datatype needed to specify type of image data
  if (! object@datatype %in% convert.datatype()) {
    retval <- c(retval, "datatype not recognized\n")
  }
  ## bitpix should correspond correctly to datatype
  if (object@bitpix != convert.bitpix()[[convert.datatype(object@datatype)]]) {
    retval <- c(retval, "bitpix does not match the datatype\n")
  }
  ## dim should be non-zero for dim[1] dimensions
  if (! all(object@"dim_"[indices] > 0)) {
    retval <- c(retval, "dim[1]/dim mismatch\n")
  }
  ## dim should be one for all >dim[1] dimensions
  if (! all(object@"dim_"[(max(indices) + 1):8] == 1)) {
    retval <- c(retval, "all dim elements > dim[1] must be 1\n")
  }
  ## number of data dimensions should match dim[1]
  if (length(indices) != length(dim(object@.Data))) {
    retval <- c(retval, "dim[1]/img mismatch\n")
  }
  ##
  if (object@"cal_min" != min(object@.Data, na.rm=TRUE) ||
      object@"cal_max" != max(object@.Data, na.rm=TRUE)) {
    retval <- c(retval, "range(img) != c(cal_min,cal_max)\n")
  }
  ## pixdim[0] is required when qform_code != 0
  if (object@"qform_code" != 0 && pixdim(object)[1] == 0) {
    retval <- c(retval, "pixdim[1] is required\n")
  }
  ## pixdim[n] required when dim[n] is required
  if (! all(object@"dim_"[indices] > 0 & pixdim(object)[indices] > 0)) {
    retval <- c(retval, "dim/pixdim mismatch\n")
  }
  ## data dimensions should match dim
  if (! isTRUE(all.equal(object@"dim_"[indices], dim(object@.Data)))) {
    retval <- c(retval, "dim/img mismatch\n")
  }
  ## vox_offset required for an "n+1" header
  if (object@"magic" == "n+1" && object@"vox_offset" == 0) {
    retval <- c(retval, "vox_offset required when magic=\"n+1\"\n")
  }
  ## magic must be "ni1\0" or "n+1\0"
  if (! (object@"magic" == "n+1" || object@"magic" == "ni1")) {
    retval <- c(retval, "magic != \"n+1\" and magic != \"ni1\"\n")
  }
  if (is.null(retval)) {
    return(TRUE)
  } else {
    return(retval)
  }
})

#############################################################################
## setClass("niftiExtension")
#############################################################################
#' @title Class "niftiExtension"
#'
#' @description An extension of the NIfTI class that allows \dQuote{extensions}
#' that conform to the NIfTI data standard.
#'
#' @name niftiExtension-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("niftiExtension", data, dim, dimnames, ...)}.
#' @author Andrew Thornton \email{zeripath@@users.sourceforge.net}
#' @seealso \code{\linkS4class{nifti}}, \code{\linkS4class{niftiAuditTrail}}
#' @references NIfTI-1\cr \url{http://nifti.nimh.nih.gov/}
#' @keywords classes
#' @examples
#' showClass("niftiExtension")
#' @export
setClass("niftiExtension",
         representation(extensions="list"),
         prototype(extensions=list()),
         contains="nifti")

#############################################################################
## setClass("niftiAuditTrail")
#############################################################################
#' @title Class "niftiAuditTrail"
#'
#' @description An extension of the NIfTI class that adds an audit trail in XML format.
#'
#' @name niftiAuditTrail-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("niftiAuditTrail", data, dim, dimnames, ...)}.
#' @author Andrew Thornton \email{zeripath@@users.sourceforge.net}
#' @seealso \code{\linkS4class{nifti}}, \code{\linkS4class{niftiExtension}}
#' @references NIfTI-1\cr \url{http://nifti.nimh.nih.gov/}
#' @section Methods:
#'   \describe{
#'     \item{show}{\code{signature(object = "niftiAuditTrail")}: prints out a
#'                 summary of the imaging data.}
#'   }
#' @keywords classes
#' @examples
#'
#' showClass("niftiAuditTrail")
#'
#' @rdname niftiAuditTrail-class
#' @export
setClass("niftiAuditTrail",
         representation(trail="ANY"),
         prototype(trail=newAuditTrail()),
         contains="niftiExtension")

#############################################################################
## setValidity("niftiExtension")
#############################################################################
setValidity("niftiExtension", function(object) {
  ## Allegedly setValidity will always check for superclasses.
  ## So we need only check that the list is empty or only contains
  ## niftiExtensionSections and check the validity of each of those
  retval <- NULL
  validSection <- getValidity(getClassDef("niftiExtensionSection"))
  lapply(object@"extensions",
         function(x) {
           if (! is(x, "niftiExtensionSection")) {
             retval <<- c(retval, paste("@extensions list contains non-niftiExtensionSection element:", class(x)))
           } else {
             if (! validSection(x)) {
               retval <<- c(retval, validSection(x))
             }
           }
         })
  if (is.null(retval)) {
    return(TRUE)
  } else {
    return(retval)
  }
})

#############################################################################
## setValidity("niftiExtensionSection")
#############################################################################
setValidity("niftiExtensionSection", function(object) {
  retval <- NULL
  if (object@esize %% 16 != 0) {
    retval <- c(retval, "esize is not a multiple of 16")
  }
  if ((object@esize - 8) < nchar(object@edata, type="bytes")) {
    retval <- c(retval, "esize is too small for the data contained within the section")
  }
  if (is.null(retval)) {
    return(TRUE)
  } else {
    return(retval)
  }
})

## setGeneric("img", function(object) { standardGeneric("img") })
## setMethod("img", "nifti", function(object) { object@.Data })
## setGeneric("img<-", function(x, value) { standardGeneric("img<-") })
## setReplaceMethod("img", signature(x="nifti", value="array"),
##           function(x, value) {
##             x@.Data <- value
##             x
##           })

## setGeneric("hdr", function(object) { standardGeneric("hdr") })
## setMethod("hdr", "nifti", # signature(object="nifti", name="ANY"),
##           function(object) { object@"descrip" })
## setGeneric("hdr<-", function(x, name, value) { standardGeneric("hdr<-") })
## setReplaceMethod("hdr", signature(x="nifti", name="character", value="ANY"),
##           function(x, name, value) {
##             x@name <- value
##             validNIfTI <- getValidity(getClassDef("nifti"))
##             validNIfTI(x)
##             x
##           })

#############################################################################
## nifti()
#############################################################################
#' @name nifti
#' @title Constructor for NIfTI
#'
#' @description Constructor for NIfTI class objects.
#'
#' @aliases nifti
#' @param img is a multidimensional array of data.
#' @param dim is the dimension of the data (default = \code{missing}).
#' @param datatype is an integer that denotes the type of data contained in
#' each voxel.  See \code{convert.datatype} or the NIfTI documentation for more
#' details.
#' @param cal.min allows user-specified minimum value in the array
#' (visualization purposes only).
#' @param cal.max allows user-specified minimum value in the array
#' (visualization purposes only).
#' @param pixdim allows user-specified pixel dimension vector (length = 8).
#' @param \dots allows for additional \sQuote{slots} to be specified.
#' @return An object of class \code{nifti}.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\linkS4class{nifti}}, \code{\link{anlz}},
#' \code{\link{convert.datatype}}
#' @references NIfTI-1\cr \url{http://nifti.nimh.nih.gov/}
#' @examples
#'
#' options("niftiAuditTrail"=FALSE)
#'
#' nim <- nifti() # default
#' nim
#' nim <- nifti(datatype=4) # 2-byte integers
#' nim
#' @export
nifti <- function(img=array(0, dim=rep(1,4)), dim, datatype=2,
                  cal.min=NULL, cal.max=NULL, pixdim=NULL, ...) {
  ## Set dimensions
  if (missing(dim)) {
    if (is.array(img)) {
      dim <- base::dim(img)
    } else {
      dim <- c(1, length(img))
    }
  }
  ld <- length(dim)
  ## Create "dim" and "pixdim" slots
  x <- rep(1, 8)
  x[1] <- length(dim)
  y <- rep(0.0, 8)
  for (i in 2:length(x)) {
    x[i] <- ifelse(is.na(dim(img)[i-1]), 1, dim(img)[i-1])
    y[i] <- ifelse(is.na(dim(img)[i-1]) || is.null(pixdim), 1.0, pixdim[i])
  }
  ## min/max values for visualization
  if (is.null(cal.max)) {
    cal.max <- as.numeric(max(img, na.rm=TRUE))
  }
  if (is.null(cal.min)) {
    cal.min <- as.numeric(min(img, na.rm=TRUE))
  }
  ## Set datatype
  switch(as.character(datatype),
         "2" = bitpix <- 8,
         "4" = bitpix <- 16,
         "8" = bitpix <- 32,
         "16" = bitpix <- 32,
         "64" = bitpix <- 64,
         "512" = bitpix <- 16,
         stop(paste("Data type", datatype, "unsupported."))
         )
  ## Create the object
  niftiClass <- "nifti"
  if (getOption("niftiAuditTrail")) {
    niftiClass <- "niftiAuditTrail"
  }
  obj <- new(niftiClass, .Data=array(img, dim=dim), "dim_"=x, "pixdim"=y,
             "cal_max"=cal.max, "cal_min"=cal.min, "datatype"=datatype,
             "bitpix"=bitpix, ...)
  if (getOption("niftiAuditTrail")) {
    audit.trail(obj) <- niftiAuditTrailCreated(call=match.call())
  }
  validNIfTI <- getValidity(getClassDef("nifti"))
  validNIfTI(obj)
  return(obj)
}

#' @title Extract or Replace NIfTI Audit Trail
#'
#' @description Operators that act on the audit trail (XML) in the NIfTI header.
#'
#' @name audit.trail-methods
#' @aliases audit.trail-methods audit.trail,nifti-method audit.trail
#' audit.trail<-,nifti-method audit.trail<-
#' @docType methods
#' @param object is of class \code{nifti}.
#' @param value Value to assign to trail slot
#' @section Methods:
#' \describe{
#' \item{object = "nifti"}{Extract or replace NIfTI audit trail.}
#' }
#' @author Andrew Thornton \email{zeripath@@users.sourceforge.net}
#' @keywords methods
#' @export
#' @rdname audit_trail-methods
setGeneric("audit.trail", function(object) { standardGeneric("audit.trail") })
#' @export
#' @rdname audit_trail-methods
#' @aliases audit.trail,nifti-method
setMethod("audit.trail", "nifti",
          function(object) {
            if (getOption("niftiAuditTrail") &&
                is(object, "niftiAuditTrail")) {
              object@"trail"
            } else {
              NULL
            }
          })
#' @export
#' @rdname audit_trail-methods
setGeneric("audit.trail<-",
           function(object, value) { standardGeneric("audit.trail<-") })
#' @export
#' @rdname audit_trail-methods
setReplaceMethod("audit.trail", "nifti",
                 function(object, value) {
                   if (getOption("niftiAuditTrail")) {
                     if (!is(object, "niftiAuditTrail")) {
                       object <- as(object, "niftiAuditTrail")
                     }
                     object@"trail" <- value
                   }
                   return(object)
                 })
#' @export
setReplaceMethod("[",
                 signature(x="nifti", i="missing", j="missing", value="array"),
                 function(x, value) {
                   x <- as.nifti(value, x)
                   validNIfTI <- getValidity(getClassDef("nifti"))
                   validNIfTI(x)
                   return(x)
                 })
#' @export
setReplaceMethod("[", signature(x="nifti", i="ANY", j="missing", value="ANY"),
                 function(x, i, value) {
                   ## For some reason this line is slow; I don't understand it
                   x@.Data[i] <- value
                   if (any(value < x@"cal_min", na.rm=TRUE)) {
                     x@"cal_min" <- min(value, na.rm=TRUE)
                   }
                   if (any(value > x@"cal_max", na.rm=TRUE)) {
                     x@"cal_max" <- max(value, na.rm=TRUE)
                   }
                   audit.trail(x) <-
                     niftiAuditTrailEvent(x, "modification", call=sys.call(-3),
                                          comment=paste("Non-numeric replace ["))
                   return(x)
                 })
#' @export
setReplaceMethod("[",
                 signature(x="nifti", i="numeric", j="missing", value="ANY"),
                 function(x, i, value) {
                   ## For some reason this line is slow; I don't understand it
                   x@.Data[i] <- value
                   if (any(value < x@"cal_min", na.rm=TRUE)) {
                     x@"cal_min" <- min(value, na.rm=TRUE)
                   }
                   if (any(value > x@"cal_max", na.rm=TRUE)) {
                     x@"cal_max" <- max(value, na.rm=TRUE)
                   }
                   audit.trail(x) <- niftiAuditTrailEvent(x, "modification",
                                                          call=sys.call(-3))
                   return(x)
                 })
#' @export
setReplaceMethod("[", signature(x="nifti", i="ANY", j="ANY", value="ANY"),
                 function(x, i, j, ..., value) {
                   ## For some reason this line is slow; I don't understand it
                   x@.Data[i,j,...] <- value
                   audit.trail(x) <-
                     niftiAuditTrailEvent(x, "modification", call=sys.call(-3),
                                          comment=paste("Non-numeric replace ["))
                   return(x)
                 })
#' @export
setReplaceMethod("[",
                 signature(x="nifti", i="numeric", j="numeric", value="ANY"),
                 function(x, i, j, ..., value) {
                   ## For some reason this line is slow; I don't understand it
                   x@.Data[i,j,...] <- value
                   if (any(value < x@"cal_min", na.rm=TRUE)) {
                     x@"cal_min" <- min(value, na.rm=TRUE)
                   }
                   if (any(value > x@"cal_max", na.rm=TRUE)) {
                     x@"cal_max" <- max(value, na.rm=TRUE)
                   }
                   audit.trail(x) <-
                     niftiAuditTrailEvent(x, "modification",
                                          comment=paste("[", paste(i, j, ..., sep=", "), "] <- ", value, sep=""))
                   return(x)
                 })


#############################################################################
## sform() accessor function to srow_*
#############################################################################
#' @title Extract NIfTI 3D Image Orientation
#'
#' @description Methods that act on the \dQuote{qform} and \dQuote{sform} information in the
#' NIfTI header.
#'
#' @name orientation-methods
#' @aliases qform-methods qform,nifti-method qform sform-methods
#' sform,nifti-method sform
#' @docType methods
#' @param object is an object of class \code{nifti}.
#' @section Methods:
#' \describe{
#' \item{object = "nifti"}{Extract or replace NIfTI description.}
#' }
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @keywords methods
#' @examples
#'
#' \dontrun{
#' url <- "http://nifti.nimh.nih.gov/nifti-1/data/avg152T1_LR_nifti.nii.gz"
#' urlfile <- file.path(system.file("nifti", package="oro.nifti"),
#'                      "mniLR.nii.gz")
#' download.file(url, urlfile, quiet=TRUE)
#' }
#' urlfile <- file.path(system.file("nifti", package="oro.nifti"),
#'                      "mniLR.nii.gz")
#' mniLR <- readNIfTI(urlfile)
#' sform(mniLR)
#' @export
#' @rdname orientation-methods
setGeneric("sform", function(object) { standardGeneric("sform") })
#' @export
#' @rdname orientation-methods
setMethod("sform", "nifti",
          function(object) {
            matrix(c(object@"srow_x", object@"srow_y", object@"srow_z"),
                   ncol=4, byrow=TRUE)
          })

#############################################################################
## qform() accessor function to quatern_*, qoffset_*
#############################################################################
#' @export
#' @rdname orientation-methods
setGeneric("qform", function(object) { standardGeneric("qform") })
#' @export
#' @rdname orientation-methods
setMethod("qform", "nifti", function(object) { quaternion2mat44(object) })
