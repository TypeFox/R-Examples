##
##
## Copyright (c) 2009-2013 Brandon Whitcher and Volker Schmid
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
##     * the names of the authors may not be used to endorse or promote
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
## $Id: convert_nifti.R 333 2010-02-02 17:43:52Z bjw34032 $
##

############################################################################
## Conversion subroutines
############################################################################
#' @name Convert NIfTI Codes
#' @title Convert NIfTI Codes
#' 
#' @description Codes that appear in the ANALYZE header are mapped to 
#' meaningful chartacter strings.
#' 
#' @details \code{switch} statements are used to map a numeric code to the 
#' appropriate string.
#' 
#' @aliases convert.bitpix convert.datatype convert.intent convert.form
#' @aliases convert.units convert.slice
#' @param bitpix is the bit-per-pixel code.
#' @param datatype.code defines data type.
#' @param intent.code is the NIfTI intent code.
#' @param form.code is the \eqn{(x,y,z)} coordinate system.
#' @param units is the units of pixdim[1..4].
#' @param inverse is a logical value that denotes the direction of unit
#' conversion.
#' @param slice.code is the slice timing order.
#' @return A character string.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @references Neuroimaging Informatics Technology Initiative (NIfTI)\cr
#' \url{http://nifti.nimh.nih.gov/}
#' @examples
#' 
#' ##  4 = SIGNED_SHORT
#' convert.datatype.anlz(4)
#' ## 16 = FLOAT
#' convert.datatype.anlz(16)
#' ##  2 = "saggital unflipped"
#' convert.orient.anlz(2)
#' ##  4 = "coronal flipped"
#' convert.orient.anlz(4)
#' 
#' @rdname convert_nifti
#' @export
convert.bitpix <- function(bitpix=NULL) {
  nifti.bitpix <- list(UINT8 = 8,
                       INT8 = 8,
                       UINT16 = 16,
                       INT16 = 16,
                       UINT32 = 32,
                       INT32 = 32,
                       UINT64 = 64,
                       INT64 = 64,
                       FLOAT32 = 32,
                       FLOAT64 = 64,
                       FLOAT128 = 128,
                       COMPLEX64 = 64,
                       COMPLEX128 = 128,
                       COMPLEX256 = 256,
                       RGB24 = 24,
                       RGBA32 = 32)
  if (is.null(bitpix)) {
    nifti.bitpix
  } else {
    names(which(nifti.bitpix == bitpix))
  }
}
#' @rdname convert_nifti
#' @export
convert.datatype <- function(datatype.code=NULL) {
  nifti.datatype <- list(UINT8 = 2,
                         INT16 = 4,
                         INT32 = 8,
                         FLOAT32 = 16,
                         COMPLEX64 = 32,
                         FLOAT64 = 64,
                         RGB24 = 128,
                         INT8 = 256,
                         UINT16 = 512,
                         UINT32 = 768,
                         INT64 = 1024,
                         UINT64 = 1280,
                         FLOAT128 = 1536,
                         COMPLEX128 = 1792,
                         COMPLEX256 = 2048,
                         RGBA32 = 2304)
  if (is.null(datatype.code)) {
    nifti.datatype
  } else {
    names(which(nifti.datatype == as.numeric(datatype.code)))
  }
}
#' @rdname convert_nifti
#' @export
convert.intent <- function(intent.code=NULL) {
  nifti.intent <- list(None = 0,
                       Correl = 2,
                       Ttest = 3,
                       Ftest = 4,
                       Zscore = 5,
                       Chisq = 6,
                       Beta = 7,
                       Binom = 8,
                       Gamma = 9,
                       Poisson = 10,
                       Normal = 11,
                       Ftest_Nonc = 12,
                       Chisq_Nonc = 13,
                       Logistic = 14,
                       Laplace = 15,
                       Uniform = 16,
                       Ttest_Nonc = 17,
                       Weibull = 18,
                       Chi = 19,
                       Invgauss = 20,
                       Extval = 21,
                       Pval = 22,
                       Logpval = 23,
                       Log10pval = 24,
                       Estimate = 1001,  # estimate of some parameter
                       Label = 1002,     # index into some set of labels
                       Neuroname = 1003, # index into the NeuroNames labels set
                       Genmatrix = 1004, # M x N matrix at each voxel
                       Symmatrix = 1005, # N x N symmetric matrix at each voxel
                       Dispvect = 1006,  # a displacement field
                       Vector = 1007,    # a displacement vector
                       Pointset = 1008,  # a spatial coordinate
                       Triangle = 1009,  # triple of indexes
                       Quaternion = 1010, # a quaternion
                       Dimless = 1011)   # Dimensionless value - no params
  if (is.null(intent.code)) {
    nifti.intent
  } else {
    names(which(nifti.intent == as.numeric(intent.code)))
  }
}
#' @rdname convert_nifti
#' @export
convert.form <- function(form.code) {
  switch(as.character(form.code),
         "0" = "Unknown",       # Arbitrary coordinates (Method 1)
         "1" = "Scanner_Anat",  # Scanner-based anatomical coordinates
         "2" = "Aligned_Anat",  # Coordinates aligned to another file's,
                                # or to anatomical "truth"
         "3" = "Talairach",     # Coordinates aligned to Talairach-
                                # Tournoux Atlas; (0,0,0)=AC, etc.
         "4" = "MNI_152")       # MNI 152 normalized coordinates
}
#' @rdname convert_nifti
#' @export
convert.units <- function(units, inverse=FALSE) {
  if (inverse) {
    switch(units,
           Unknown = 0,           # unspecified units
           meter = 1,             # meters
           mm = 2,                # millimeters
           micron = 3,            # micrometers
           sec = 8,               # seconds
           msec = 16,             # milliseconds
           usec = 24,             # microseconds
           Hz = 32,               # Hertz
           ppm = 40,              # parts per million
           rads = 48)             # radians per second
  } else {
    switch(as.character(units),
           "0" = "Unknown",       # unspecified units
           "1" = "meter",         # meters
           "2" = "mm",            # millimeters
           "3" = "micron",        # micrometers
           "8" = "sec",           # seconds
           "16" = "msec",         # milliseconds
           "24" = "usec",         # microseconds
           "32" = "Hz",           # Hertz
           "40" = "ppm",          # parts per million
           "48" = "rads")         # radians per second
  }
}
#' @rdname convert_nifti
#' @export
convert.slice <- function(slice.code) {
  switch(as.character(slice.code),
         "0" = "Unknown",
         "1" = "Seq_Inc",       # sequential increasing
         "2" = "Seq_Dec",       # sequential decreasing
         "3" = "Alt_Inc",       # alternating increasing
         "4" = "Alt_Dec",       # alternating decreasing
         "5" = "Alt_Inc2",      # alternating increasing #2
         "6" = "Alt_Dec2")      # alternating decreasing #2
}

############################################################################
## Bitwise conversion subroutines
############################################################################
#' Bitwise Conversion Subroutines
#' 
#' Units of spatial and temporal dimensions, and MRI-specific spatial and
#' temporal information.
#' 
#' The functions \code{xyzt2space} and \code{xyzt2time} can be used to mask off
#' the undesired bits from the \kbd{xyzt_units} fields, leaving \dQuote{pure}
#' space and time codes.\cr
#' \url{http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/xyzt_units.html}
#' 
#' The functions \code{dim2freq}, \code{dim2phase}, and \code{dim2slice} can be
#' used to extract values from the \kbd{dim_info} byte.\cr
#' \url{http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/dim_info.html}
#' 
#' @aliases xyzt2space xyzt2time space.time2xyzt dim2freq dim2phase dim2slice
#' @param xyzt represents the units of pixdim[1..4] in the NIfTI header.
#' @param ss is the character string of spatial units.  Valid strings are:
#' \dQuote{Unknown}, \dQuote{meter}, \dQuote{mm} and \dQuote{micron}.
#' @param tt is the character string of temporal units.  Valid strings are:
#' \dQuote{sec}, \dQuote{msec}, \dQuote{usec}, \dQuote{Hz}, \dQuote{ppm} and
#' \dQuote{rads}.
#' @param di represents MRI slice ordering in the NIfTI header.
#' @return For \kbd{diminfo}: the frequency, phase and slice dimensions encode
#' which spatial dimension (1,2, or 3) corresponds to which acquisition
#' dimension for MRI data.  For \kbd{xyzt_units}: the codes are used to
#' indicate the units of pixdim.  Dimensions 1,2,3 are for x,y,z; dimension 4
#' is for time (t).
#' @author B. Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{convert.units}}, \code{\link{convert.slice}}
#' @references Neuroimaging Informatics Technology Initiative (NIfTI)\cr
#' \url{http://nifti.nimh.nih.gov/}
#' @keywords misc
#' @rdname bitwise
#' @import bitops
#' @export
xyzt2space <- function(xyzt) {
  ## define XYZT_TO_SPACE(xyzt) ( (xyzt) & 0x07 )
  bitops::bitAnd(xyzt, 7)
}
#' @rdname bitwise
#' @export
xyzt2time <- function(xyzt) {
  ## define XYZT_TO_TIME(xyzt) ( (xyzt) & 0x38 )
  bitops::bitAnd(xyzt, 56)
}
#' @rdname bitwise
#' @export
space.time2xyzt <- function(ss, tt) {
  ## define SPACE_TIME_TO_XYZT(ss,tt) (  (((char)(ss)) & 0x07)   \
  ##                                   | (((char)(tt)) & 0x38) )
  ss <- bitops::bitAnd(convert.units(ss, inverse=TRUE), 7)
  tt <- bitops::bitAnd(convert.units(tt, inverse=TRUE), 56)
  ss + tt
}
#' @rdname bitwise
#' @export
dim2freq <- function(di) {
  ## define DIM_INFO_TO_FREQ_DIM(di) ( ((di)     ) & 0x03 )
  bitops::bitAnd(di, 3)
}
#' @rdname bitwise
#' @export
dim2phase <- function(di) {
  ## define DIM_INFO_TO_PHASE_DIM(di) ( ((di) >> 2) & 0x03 )
  bitops::bitAnd(bitops::bitShiftR(di, 2), 3)
}
#' @rdname bitwise
#' @export
dim2slice <- function(di) {
  ## define DIM_INFO_TO_SLICE_DIM(di) ( ((di) >> 4) & 0x03 )
  bitops::bitAnd(bitops::bitShiftR(di, 4), 3)
}

############################################################################
## as.nifti()
############################################################################
#' @name as.nifti
#' @title as.nifti
#' 
#' @description Internal function that converts multidimensional arrays to 
#' NIfTI class objects.
#' 
#' @aliases as.nifti
#' @param from is the object to be converted.
#' @param value is the \code{anlz} class object to use as a template for
#' various NIfTI header information.
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @return An object of class \code{nifti}.
#' @author Andrew Thornton \email{zeripath@@users.sourceforge.net}.\cr
#' Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @export
as.nifti <- function(from, value=NULL, verbose=FALSE) {
  anlz.as.nifti <- function(from, value=nifti()) {
    ## So what kind of thing do we keep?
    slots <- c("dim_", "datatype", "bitpix", "pixdim", "descrip",
               "aux_file", ".Data")
    for (s in 1:length(slots)) {
      slot(value, slots[s]) <- slot(from, slots[s])
    }
    ## sapply(slots, function(x) { slot(value, x) <<- slot(from, x); NULL })
    calset <- !(from@"cal_max" == 0 && from@"cal_min" == 0)
    slot(value, "cal_max") <- ifelse(calset, from@"cal_min", from@"glmax")
    slot(value, "cal_min") <- ifelse(calset, from@"cal_max", from@"glmin")
    return(value)
    ## The below code should apply the orient code as per 
    ## http://eeg.sourceforge.net/ANALYZE75.pdf and 
    ## http://eeg.sourceforge.net/AnalyzeDirect_ANALYZEformat.pdf
    ## However the NIfTI website says that this field is not often
    ##   set properly so I am unsure whether to apply this transform
    ##
    ##R<-diag(1, 4)
    ##R[,1]<- -R[,1] # as i is by default Left in the ANALYZE standard
    ##switch(from@"orient",
    ## orient:   slice orientation for this dataset.
    ##      0         transverse unflipped (i,j,k)=(L,A,S)
    #0=R,
    ##      1         coronal unflipped    (i,j,k)=(L,S,A)
    #1=R[,2:3] <- R[,3:2],
    ##      2         sagittal unflipped   (i,j,k)=(A,S,L)
    #2=R[,c(1,3)] <- R[,c(3,1)],
    ##      3         transverse flipped   (i,j,k)=(R,A,S)
    #3=R[,1]<--R[,1],
    ##      4         coronal flipped	    (i,j,k)=(R,S,A)
    #4=R[,1]<--R[,1];R[,2:3] <- R[,3:2],
    ##      5         sagittal flipped	    (i,j,k)=(A,S,R)
    #5=R[,1]<--R[,1];R[,c(1,3)] <- R[,c(3,1)],
    #stop("Unknown orientation: ",from@orient))
    #value@"xyzt_units" <- from@"vox_units"
    #value@"sform_code" <- 1
    #value@"srow_x" <- R[1,]
    #value@"srow_y" <- R[2,]
    #value@"srow_z" <- R[3,]
  }
  integertype <- function(from) {
    integer.ranges <- list("UINT8" = c(0, 2^8-1),
                           "INT16" = c(-2^15, 2^15-1),
                           "INT32" = c(-2^31, 2^31-1),
                           "INT8" = c(-2^7, 2^7-1),
                           "INT64" = c(-2^63, 2^63-1),
                           "UINT64" = c(0, 2^64))
    fromRange <- range(from)
    for (i in 1:length(integer.ranges)) {
      if (fromRange[1] >= integer.ranges[[i]][1] &&
	  fromRange[2] <= integer.ranges[[i]][2]) {
	return(names(integer.ranges)[i])
      }
    }
    warning("Range too large to be kept as integer, forcing float")
    floattype(from)
  }
  floattype <- function(from) {
    return("FLOAT32")
  }
  if (is.null(value)) {
    nim <- nifti()
  } else {
    nim <- value
  }
  if (is(from, "anlz")) {
    nim <- anlz.as.nifti(from)
  } else {
    if (is.array(from)) {
      ## Determine a sensible datatype
      dataClass <- class(from[1])
      datatypeString <- switch(dataClass,
                               logical = integertype(from),
                               integer = integertype(from),
                               numeric = floattype(from),
                               stop("Can't transform data in from: ",
                                    class(from[1])))
      nim@"data_type" <- datatypeString
      nim@"datatype" <- convert.datatype()[[datatypeString]]
      nim@"bitpix" <- convert.bitpix()[[datatypeString]]
      nim@"cal_min" <- min(from, na.rm=TRUE)
      nim@"cal_max" <- max(from, na.rm=TRUE)
      nim@"dim_" <- c(length(dim(from)), dim(from))
      if (length(nim@"dim_") < 8) {
        nim@"dim_" <- c(nim@"dim_", rep(1, 8 - length(nim@"dim_")))
      }
      nim@.Data <- from
      if (getOption("niftiAuditTrail") && is(nim, "niftiAuditTrail")) {
        audit.trail(nim) <-
          niftiAuditTrailCreated(history=nim, call=match.call())
      }
    } else {
      if (is.list(from)) {
        nim <- lapply(from, function(x) as.nifti(x, value))
        lapply(names(from),
               function(x) {
                 if (is.nifti(nim[[x]])) {
                   nim[[x]]@"intent_code" <<- convert.intent()[["Estimate"]]
                   nim[[x]]@"intent_name" <<- substr(x, 1, 15)
                   audit.trail(nim[[x]]) <<-
                     niftiAuditTrailEvent(nim[[x]],
                                          type="Intent Changed",
                                          comment=paste("Parameter Estimate:", x))
                 }
               })
      } else {
        if (verbose) {
          warning("Cannot convert class =", class(from), "to nifti object")
        }
        nim <- from
      }
    }
  }
  return(nim)
}
