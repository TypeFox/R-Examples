##
##
## Copyright (c) 2005-2010 Karsten Tabelow (WIAS)
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
## $Id: afniS4.R 2010-12-17 15:17:07 $
##

#############################################################################
## setClass("afni")
#############################################################################
#' Class "afni"
#' 
#' The AFNI class for medical imaging data.
#' 
#' @name afni-class
#' @aliases afni-class show,afni-method
#' @param object An object of class \code{afni}.
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("afni", data, dim, dimnames, ...)}.
#' @author Karsten Tabelow \email{karsten.tabelow@@wias-berlin.de}
#' @seealso \code{\linkS4class{nifti}}, \code{\linkS4class{anlz}}
#' @references AFNI\cr
#' \url{http://afni.nimh.nih.gov/pub/dist/src/README.attributes}
#' @keywords classes
#' @examples
#' 
#' showClass("afni")
#' @section Slots:
#'   \describe{
#'     \item{\code{.Data}:}{Object of class \code{"array"} contains the imaging data}
#'     \item{\code{DATASET_RANK}:}{Object of class \code{"integer"}}
#'     \item{\code{DATASET_DIMENSIONS}:}{Object of class \code{"integer"}}
#'     \item{\code{TYPESTRING}:}{Object of class \code{"character"}}
#'     \item{\code{SCENE_DATA}:}{Object of class \code{"integer"}}
#'     \item{\code{ORIENT_SPECIFIC}:}{Object of class \code{"integer"}}
#'     \item{\code{ORIGIN}:}{Object of class \code{"numeric"}}
#'     \item{\code{DELTA}:}{Object of class \code{"numeric"}}
#'     \item{\code{TAXIS_NUMS}:}{Object of class \code{"integer"}}
#'     \item{\code{TAXIS_FLOATS}:}{Object of class \code{"numeric"}}
#'     \item{\code{TAXIS_OFFSETS}:}{Object of class \code{"numeric"}}
#'     \item{\code{IDCODE_STRING}:}{Object of class \code{"character"}}
#'     \item{\code{IDCODE_DATE}:}{Object of class \code{"character"}}
#'     \item{\code{BYTEORDER_STRING}:}{Object of class \code{"character"}}
#'     \item{\code{BRICK_STATS}:}{Object of class \code{"numeric"}}
#'     \item{\code{BRICK_TYPES}:}{Object of class \code{"integer"}}
#'     \item{\code{BRICK_FLOAT_FACS}:}{Object of class \code{"numeric"}}
#'     \item{\code{BRICK_LABS}:}{Object of class \code{"character"}}
#'     \item{\code{BRICK_STATAUX}:}{Object of class \code{"numeric"}}
#'     \item{\code{STAT_AUX}:}{Object of class \code{"numeric"}}
#'     \item{\code{HISTORY_NOTE}:}{Object of class \code{"character"}}
#'     \item{\code{NOTES_COUNT}:}{Object of class \code{"integer"}}
#'     \item{\code{NOTE_NUMBER}:}{Object of class \code{"character"}}
#'     \item{\code{TAGALIGN_MATVEC}:}{Object of class \code{"numeric"}}
#'     \item{\code{VOLREG_MATVEC}:}{Object of class \code{"array"}}
#'     \item{\code{VOLREG_ROTCOM}:}{Object of class \code{"character"}}
#'     \item{\code{VOLREG_CENTER_OLD}:}{Object of class \code{"numeric"}}
#'     \item{\code{VOLREG_CENTER_BASE}:}{Object of class \code{"numeric"}}
#'     \item{\code{VOLREG_ROTPARENT_IDCODE}:}{Object of class \code{"character"}}
#'     \item{\code{VOLREG_ROTPARENT_NAME}:}{Object of class \code{"character"}}
#'     \item{\code{VOLREG_GRIDPARENT_IDCODE}:}{Object of class \code{"character"}}
#'     \item{\code{VOLREG_GRIDPARENT_NAME}:}{Object of class \code{"character"}}
#'     \item{\code{VOLREG_INPUT_IDCODE}:}{Object of class \code{"character"}}
#'     \item{\code{VOLREG_INPUT_NAME}:}{Object of class \code{"character"}}
#'     \item{\code{VOLREG_BASE_IDCODE}:}{Object of class \code{"character"}}
#'     \item{\code{VOLREG_BASE_NAME}:}{Object of class \code{"character"}}
#'     \item{\code{VOLREG_ROTCOM_NUM}:}{Object of class \code{"integer"}}
#'     \item{\code{IDCODE_ANAT_PARENT}:}{Object of class \code{"character"}}
#'     \item{\code{TO3D_ZPAD}:}{Object of class \code{"integer"}}
#'     \item{\code{IDCODE_WARP_PARENT}:}{Object of class \code{"character"}}
#'     \item{\code{WARP_TYPE}:}{Object of class \code{"integer"}}
#'     \item{\code{WARP_DATA}:}{Object of class \code{"numeric"}}
#'     \item{\code{MARKS_XYZ}:}{Object of class \code{"numeric"}}
#'     \item{\code{MARKS_LAB}:}{Object of class \code{"character"}}
#'     \item{\code{MARKS_HELP}:}{Object of class \code{"character"}}
#'     \item{\code{MARKS_FLAGS}:}{Object of class \code{"integer"}}
#'     \item{\code{TAGSET_NUM}:}{Object of class \code{"integer"}}
#'     \item{\code{TAGSET_FLOATS}:}{Object of class \code{"numeric"}}
#'     \item{\code{TAGSET_LABELS}:}{Object of class \code{"character"}}
#'     \item{\code{LABEL_1}:}{Object of class \code{"character"}}
#'     \item{\code{LABEL_2}:}{Object of class \code{"character"}}
#'     \item{\code{DATASET_NAME}:}{Object of class \code{"character"}}
#'     \item{\code{DATASET_KEYWORDS}:}{Object of class \code{"character"}}
#'     \item{\code{BRICK_KEYWORDS}:}{Object of class \code{"character"}}
#'   }
#' 
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
#'   @export
#'   @rdname afni-class
setClass("afni", 
         representation(
           ## Mandatory attributes                 # number in R index count
           DATASET_RANK = "integer",               # 2
           DATASET_DIMENSIONS = "integer",         # 3
           TYPESTRING = "character",               # 1
           SCENE_DATA = "integer",                 # 3
           ORIENT_SPECIFIC = "integer",            # 3
           ORIGIN = "numeric",                     # 3
           DELTA = "numeric",                      # 3
           ## Mandatory if 3D+t dataset
           TAXIS_NUMS = "integer",                 # 3
           TAXIS_FLOATS = "numeric",               # 5
           TAXIS_OFFSETS = "numeric",              # TAXIS_NUMS[2]
           ## Almost mandatory attributes
           IDCODE_STRING = "character",            # 1
           IDCODE_DATE = "character",              # 1
           BYTEORDER_STRING = "character",         # 1
           BRICK_STATS = "numeric",                # 2*DATASET_RANK[2]
           BRICK_TYPES = "integer",                # DATASET_RANK[2]
           BRICK_FLOAT_FACS = "numeric",           # DATASET_RANK[2]
           BRICK_LABS = "character",               # DATASET_RANK[2]
           BRICK_STATAUX = "numeric",              # ???          
           STAT_AUX = "numeric",                   # ???
           ## Notes Attributes
           HISTORY_NOTE = "character",             # 1
           NOTES_COUNT = "integer",                # 1
           NOTE_NUMBER = "character",              # NOTES_COUNT
           ## Registration Attributes
           TAGALIGN_MATVEC = "numeric",            # 12
           VOLREG_MATVEC = "array",                # 12xDATASET_RANK[2]
           VOLREG_ROTCOM = "character",            # DATASET_RANK[2]
           VOLREG_CENTER_OLD = "numeric",          # 3?
           VOLREG_CENTER_BASE = "numeric",         # 3?
           VOLREG_ROTPARENT_IDCODE = "character",  # 1
           VOLREG_ROTPARENT_NAME = "character",    # 1
           VOLREG_GRIDPARENT_IDCODE = "character", # 1
           VOLREG_GRIDPARENT_NAME = "character",   # 1
           VOLREG_INPUT_IDCODE = "character",      # 1
           VOLREG_INPUT_NAME = "character",        # 1
           VOLREG_BASE_IDCODE = "character",       # 1
           VOLREG_BASE_NAME = "character",         # 1
           VOLREG_ROTCOM_NUM = "integer",          # 1
           ## Miscellaneous Attributes
           IDCODE_ANAT_PARENT = "character",       # 1
           TO3D_ZPAD = "integer",                  # 3
           ## Warping Attributes
           IDCODE_WARP_PARENT = "character",       # 1
           WARP_TYPE = "integer",                  # 2
           WARP_DATA = "numeric",                  # 30 (WARP_TYPE[1] == 0) or 360 (WARP_TYPE[1] == 1)
           ## Talairach Markers Attributes
           MARKS_XYZ = "numeric",                  # 3x10
           MARKS_LAB = "character",                # 10x20
           MARKS_HELP = "character",               # 10x256
           MARKS_FLAGS = "integer",                # 2
           ## Attributes for User-Defined Tags
           TAGSET_NUM = "integer",                 # 2 (TAGSET_NUM[2] = 5)
           TAGSET_FLOATS = "numeric",              # TAGSET_NUM[1]*TAGSET_NUM[2]
           TAGSET_LABELS = "character",            # TAGSET_NUM[1]
           ## Nearly Useless Attributes
           LABEL_1 = "character",                  # 1
           LABEL_2 = "character",                  # 1
           DATASET_NAME = "character",             # 1
           DATASET_KEYWORDS = "character",         # 1
           BRICK_KEYWORDS = "character"            # 1
         ),
         prototype(
           ## Mandatory attributes 
           DATASET_RANK = integer(2),
           DATASET_DIMENSIONS = integer(3),
           TYPESTRING = "",
           SCENE_DATA = integer(3),
           ORIENT_SPECIFIC = integer(3),
           ORIGIN = numeric(3),
           DELTA = numeric(3),
           ## Mandatory if 3D+t dataset
           TAXIS_NUMS = integer(3),
           TAXIS_FLOATS = numeric(5),
           TAXIS_OFFSETS = numeric()
         ),
         contains="array")

#############################################################################
## setMethod("show", "afni")
#############################################################################
#' @rdname afni-class
#' @aliases show,afni-method
setMethod("show", "afni", function(object) {
  cat("AFNI format", fill=TRUE)
  cat("  Type            :", class(object), fill=TRUE)
  cat("  Typestring      :", object@TYPESTRING, fill=TRUE)
  cat("  Scene Data      :", object@SCENE_DATA[1:3],
      " (", convert.scene(object@SCENE_DATA[2], object@TYPESTRING), ")", 
      fill=TRUE)
  cat("  Dimension       :",
      paste(c(object@DATASET_DIMENSIONS[1:3], 
              object@DATASET_RANK[2]), collapse=" x "),
      fill=TRUE)
  cat("  Pixel Dimension :",
      paste(round(object@DELTA[1:3],2), collapse=" x "),
      fill=TRUE)
})

#############################################################################
## setValidity("afni")
#############################################################################
setValidity("afni", function(object) {
  retval <- NULL
  ## test existence of mandatory attributes
  if (length(object@DATASET_RANK) < 2)
    retval <- c(retval, "missing or incomplete attribute DATASET_RANK")
  if (length(object@DATASET_DIMENSIONS) < 3)
    retval <- c(retval, "missing or incomplete attribute DATASET_DIMENSIONS")
  if (length(object@TYPESTRING) < 1)
    retval <- c(retval, "missing or incomplete attribute TYPESTRING")
  if (length(object@SCENE_DATA) < 3)
    retval <- c(retval, "missing or incomplete attribute SCENE_DATA")
  if (length(object@ORIENT_SPECIFIC) < 3)
    retval <- c(retval, "missing or incomplete attribute ORIENT_SPECIFIC")
  if (length(object@ORIGIN) < 3)
    retval <- c(retval, "missing or incomplete attribute ORIGIN")
  if (length(object@DELTA) < 3)
    retval <- c(retval, "missing or incomplete attribute DELTA")
  if (!is.null(retval)) {
    return(retval)
  }
  ## consistency check for mandatory attributes
  if (object@DATASET_RANK[1] != 3)
    retval <- c(retval, "DATASET_RANK[1] is not 3")
  if (object@DATASET_RANK[2] != dim(object@.Data)[4])
    retval <- c(retval, "DATASET_RANK[2] does not equal number of sub-bricks in .Data")
  if (any(object@DATASET_DIMENSIONS[1:3] != dim(object@.Data)[1:3]))
    retval <- c(retval, "DATASET_RANK[2] does not equal number of sub-bricks in .Data")
  if (!(object@TYPESTRING %in% c("3DIM_HEAD_ANAT",
                                 "3DIM_HEAD_FUNC",
                                 "3DIM_GEN_ANAT",
                                 "3DIM_GEN_ANAT")))
    retval <- c(retval, paste("unknown TYPESTRING", object@TYPESTRING))
  if ((object@SCENE_DATA[1] < 0) || (object@SCENE_DATA[1] > 2))
    retval <- c(retval, "view type in SCENE_DATA out of range 0, 1, 2")
  if ((object@SCENE_DATA[2] < 0) || (object@SCENE_DATA[2] > 11))
    retval <- c(retval, "func type in SCENE_DATA out of range 0, 1, ..., 11")
  if (object@SCENE_DATA[3] + 1 != which(object@TYPESTRING == c("3DIM_HEAD_ANAT", "3DIM_HEAD_FUNC", "3DIM_GEN_ANAT", "3DIM_GEN_ANAT")))
    retval <- c(retval, "malformed SCENE_DATA[3] does not match TYPESTRING")
  if (any(object@ORIENT_SPECIFIC < 0) || any(object@ORIENT_SPECIFIC > 5))
    retval <- c(retval, "ORIENT_SPECIFIC out of range 0, 1, ..., 5")
  ## more consistency checks, but certainly not all!
  ## see http://afni.nimh.nih.gov/pub/dist/src/README.attributes
  if (length(object@TAXIS_NUMS) > 0)
    if (object@TAXIS_NUMS[1] != object@DATASET_RANK[2])
      retval <- c(retval, "TAXIS_NUMS[1] does not equal DATASET_RANK[2]")
  if (length(object@TAXIS_NUMS) > 0)
    if (object@TAXIS_NUMS[2] > 0)
      if (object@TAXIS_NUMS[2] != length(object@TAXIS_OFFSETS))
        retval <- c(retval, "TAXIS_NUMS[2] does not equal length of TAXIS_OFFSETS")
  if (length(object@BYTEORDER_STRING) > 0)
    if (!(object@BYTEORDER_STRING %in% c("LSB_FIRST", "MSB_FIRST")))
      retval <- c(retval, paste("unknown BYTEORDER_STRING", object@BYTEORDER_STRING))
  if (length(object@BRICK_TYPES) > 0)
    if (object@DATASET_RANK[2] != length(object@BRICK_TYPES))
      retval <- c(retval, "DATASET_RANK[2] does not equal length of BRICK_TYPES")
  if (is.null(retval)) {
    return(TRUE)
  } else {
    return(retval)
  }
})

#############################################################################
## readAFNI()
#############################################################################
#' @title readAFNI
#' 
#' @description These functions read in the header information and
#' multidimensional array from a binary file in AFNI format into a
#' \code{\linkS4class{afni}}-class object.
#' 
#' @details The \code{readAFNI} function utilizes internal methods \code{readBin} and
#' \code{readLines} to efficiently extract information from the header and
#' binary file(s).  Compression is allowed on the BRIK file using gzip.
#' 
#' Current acceptable data types include \describe{ \item{list("INT16")}{DT
#' SIGNED SHORT (16 bits per voxel)} \item{list("FLOAT32")}{DT FLOAT (32 bits
#' per voxel)} \item{list("COMPLEX128")}{DT COMPLEX (128 bits per voxel)} }
#' 
#' @param fname is the file name of the AFNI file.
#' @param vol vector of brick numbers to be read from file.
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @param warn is a number to regularegulatete the display of warnings (default
#' = -1).  See \code{options} for more details.
#' @param call keeps track of the current function call for use in the AFNI
#' extension.
#' @return object of class \code{\linkS4class{afni}}
#' @author Karsten Tabelow \email{karsten.tabelow@@wias-berlin.de}
#' @seealso \code{\link{readANALYZE}}, \code{\link{readNIfTI}}
#' @references AFNI\cr
#' \url{http://afni.nimh.nih.gov/pub/dist/src/README.attributes}
#' @keywords file methods
#' @export
#' @examples
#' \dontrun{
#' ## Taken from the AFNI Matlab Library
#' ## http://afni.nimh.nih.gov/pub/dist/data/afni_matlab_data.tgz
#' afni.path <- system.file("afni", package="oro.nifti")
#' orig <- readAFNI(file.path(afni.path, "ARzs_CW_avvr.DEL+orig"))
#' image(orig, zlim=c(0.5,256), oma=rep(2,4))
#' orthographic(orig, zlim=c(0.5,256), oma=rep(2,4))
#' ## Taken from the AFNI installation
#' TT <- readAFNI(file.path(afni.path, "TT_N27_EZ_LR+tlrc"))
#' image(TT, zlim=c(0.5,256), oma=rep(2,4))
#' orthographic(TT, zlim=c(0.5,256), oma=rep(2,4))
#' }
#' @rdname read_afni
#' @name readAFNI
readAFNI <- function(fname, vol=NULL, verbose=FALSE, warn=-1, call=NULL) {
  if (is.null(call)) {
    call <- match.call()
  }
  ## Warnings?
  oldwarn <- getOption("warn")
  options(warn=warn)
  if (verbose) {
    cat(paste("  fname =", fname), fill=TRUE)
  }
  ## Strip any extensions
  fname <- sub("\\.gz$", "", fname)
  fname <- sub("\\.HEAD$", "", fname)
  fname <- sub("\\.head$", "", fname)
  fname <- sub("\\.BRIK$", "", fname)
  fname <- sub("\\.brik$", "", fname)

  if (! (file.exists(paste(fname, "HEAD", sep=".")) &&
         file.exists(paste(fname, "BRIK", sep="."))) &&
      ! (file.exists(paste(fname, "HEAD", sep=".")) &&
         file.exists(paste(fname, "BRIK.gz", sep=".")))) {
    stop("Header or image file(s) not found! (expect Extension HEAD/BRIK)")
  }
  ## if exist, upload
  #  if (verbose) {
  #    cat(paste("  files = ", fname, ".HEAD/BRIK", sep=""), fill=TRUE)
  #  }
  ## If uncompressed files exist, then upload!
  if ((file.exists(paste(fname, "head", sep=".")) &&
       file.exists(paste(fname, "brik", sep="."))) || 
      (file.exists(paste(fname, "HEAD", sep=".")) &&
       file.exists(paste(fname, "BRIK", sep=".")))) {      
    if (verbose) {
      cat(paste("  files = ", fname, ".HEAD/BRIK", sep=""), fill=TRUE)
    }
    aim <- .read.afni.content(fname, vol=vol, gzipped=FALSE, verbose=verbose, 
                              warn=warn, call=call)
    #    options(warn=oldwarn)
    #    return(aim)
  }
  ## If compressed files exist, then upload!
  if ((file.exists(paste(fname, "head", sep=".")) &&
       file.exists(paste(fname, "brik.gz", sep="."))) ||
      (file.exists(paste(fname, "HEAD", sep=".")) &&
       file.exists(paste(fname, "BRIK.gz", sep=".")))) {      
    if (verbose) {
      cat(paste("  files = ", fname, ".HEAD/BRIK.gz", sep=""), fill=TRUE)
    }
    aim <- .read.afni.content(fname, vol=vol, gzipped=TRUE, verbose=verbose, 
                              warn=warn, call=call)
    #    options(warn=oldwarn)
    #    return(aim)
  }
  
  #  nim <- .read.afni.content(fname, vol=vol, verbose=verbose, warn=warn,
  #                            call=call)
  options(warn=oldwarn)
  #  invisible(nim)
  invisible(aim)
}

############################################################################
############################################################################
############################################################################

.read.afni.content <- function(fname, vol, gzipped=FALSE, verbose=FALSE, 
                               warn=-1, call=NULL) {
  ## we know these files exist
  fname.head <- paste(fname, "HEAD", sep=".")
  if (gzipped) {
      fname.brik <- paste(fname, "BRIK", "gz", sep=".")
  } else {
      fname.brik <- paste(fname, "BRIK", sep=".")
  }
  ## Warnings?
  oldwarn <- getOption("warn")
  options(warn=warn)
  ## read header file
  ## Open header file
  conhead <- file(fname.head, "r")
  header <- readLines(conhead)
  close(conhead)
  types <- NULL
  args <- NULL
  counts <- NULL
  values <- NULL
  for (i in 1:length(header)) {
    if (regexpr("^type *= *", header[i]) != -1) {
      tmptype <- strsplit(header[i], " *= *")[[1]][2]
      types <- c(types, tmptype)
      args <- c(args, strsplit(header[i+1], " *= *")[[1]][2])
      tmpcounts <- as.numeric(strsplit(header[i+2], " *= *")[[1]][2])
      counts <- c(counts, tmpcounts)
      i <- i + 3
      tmpvalue <- ""
      while ((regexpr("^$", header[i]) == -1) && (i <= length(header))) {
        tmpvalue <- paste(tmpvalue,header[i])
        i <- i + 1
      }
      tmpvalue <- sub("^ +", "", tmpvalue)
      if (tmptype == "integer-attribute") {
        tmpvalue <- as.integer(strsplit(tmpvalue, " +")[[1]])
      } else if (tmptype == "float-attribute") {
        tmpvalue <- as.numeric(strsplit(tmpvalue, " +")[[1]])
      } else {
        tmpvalue <- sub("~$", "", sub("^\'", "", tmpvalue))
      }
      values <- c(values, list(value = tmpvalue))
    }        
  }
  names(values) <- args
  values$NOTE_NUMBER <- character()
  if (any(ind <- regexpr("NOTE_NUMBER_[0-9]", args) != -1)) {
    values$NOTE_NUMBER <- character(sum(ind))
    for (i in 1:length(args)) {
      arg <- args[i]
      if (length(tmp <- strsplit(arg, "NOTE_NUMBER_")[[1]]) == 2) {
        values$VOLREG_ROTCOM[as.integer(tmp[2])] <- values[i]
      }
    }
  } else {
    values$NOTE_NUMBER <- NULL
  }
  values$VOLREG_ROTCOM <- character()
  if (any(ind <- regexpr("VOLREG_ROTCOM_[0-9]", args) != -1)) {
    values$VOLREG_ROTCOM <- character(sum(ind))
    for (i in 1:length(args)) {
      arg <- args[i]
      if (length(tmp <- strsplit(arg, "VOLREG_ROTCOM_")[[1]]) == 2) {
        values$VOLREG_ROTCOM[as.integer(tmp[2]) + 1] <- values[i][[1]]
      }
    }
  } else {
    values$VOLREG_ROTCOM <- NULL
  }
  values$VOLREG_MATVEC <- numeric()
  if (any(ind <- regexpr("VOLREG_MATVEC_[0-9]", args) != -1)) {
    values$VOLREG_MATVEC <- numeric(12*sum(ind))
    for (i in 1:length(args)) {
      arg <- args[i]
      if (length(tmp <- strsplit(arg, "VOLREG_MATVEC_")[[1]]) == 2) {
        values$VOLREG_MATVEC[(12*as.integer(tmp[2])+1):(12*(as.integer(tmp[2])+1))] <- values[i][[1]]
      }
    }
    dim(values$VOLREG_MATVEC) <- c(12, sum(ind))
  } else {
    values$VOLREG_MATVEC <- NULL
  }
  ## read image data
  ddim <- values$DATASET_DIMENSIONS[1:3]
  dt <- values$DATASET_RANK[2]
  scale <- values$BRICK_FLOAT_FACS
  ## size <- file.info(fname.brik)$size / (prod(ddim) * dt)
  if (is.null(values$BRICK_TYPES[1])) {
      stop("BRICK_TYPES is unknown!")
#    what <- switch(as.character(size),
#                   "2" = "integer",
#                   "4" = "double",
#                   "16" = "complex",
#                   "integer")
  } else {
    what <- switch(as.character(values$BRICK_TYPES[1]),
                   "1" = "integer",
                   "3" = "double",
                   "5" = "complex",
                   "integer")
    size <- switch(what, 
                   "integer" = 2,
                   "double" = 4,
                   "complex" = 16)
  }
  endian <- if (regexpr("MSB", values$BYTEORDER_STRING[1]) != -1) "big" else "little"
  if (verbose) {
    cat(paste("  endianess =", endian), fill=TRUE)
  }
  if (is.null(vol)) {
    vol <- 1:dt
  }
  dataCube <- numeric(prod(ddim) * length(vol))
  if ((as.integer(size) == size) && (length(vol) > 0)) {
    kk <- 1
    if (gzipped) {
        conbrik <- gzfile(fname.brik, "rb")
    } else {
        conbrik <- file(fname.brik, "rb")
    }
#    conbrik <- file(fname.brik, "rb")
    for (k in 1:dt) {
      if (k %in% vol) {
        if (scale[k] != 0) {
          dataCube[((kk-1) * prod(ddim) + 1):(kk * prod(ddim))] <-
            readBin(conbrik,
                    what,
                    n = prod(ddim),
                    size = size,
                    signed = TRUE,
                    endian = endian) * scale[k]
        } else {
          dataCube[((kk-1) * prod(ddim) + 1):(kk * prod(ddim))] <-
            readBin(conbrik,
                    what,
                    n = prod(ddim),
                    size = size,
                    signed = TRUE,
                    endian = endian)
        }
        kk <- kk + 1
      } else {
        seek(conbrik, where=prod(ddim)*size, origin="current")
      }
    }
    close(conbrik)
  } else {
    stop("Error reading file", fname.brik,
         ": Could not detect size per voxel\n")
  }
  ## Constructor S4 object
  dim(dataCube) <- c(ddim, length(vol))
  values$DATASET_RANK[2] <- length(vol)
  if (!is.null(values$TAXIS_NUMS)) values$TAXIS_NUMS[1] <- length(vol)
  if (!is.null(values$BRICK_TYPES)) values$BRICK_TYPES <- values$BRICK_TYPES[vol]
  nim <- new("afni",
             dataCube,
             DATASET_RANK = values$DATASET_RANK,
             DATASET_DIMENSIONS = values$DATASET_DIMENSIONS,
             TYPESTRING = values$TYPESTRING,
             SCENE_DATA = values$SCENE_DATA,
             ORIENT_SPECIFIC = values$ORIENT_SPECIFIC,
             ORIGIN = values$ORIGIN,
             DELTA = values$DELTA,
             TAXIS_NUMS = if (is.null(values$TAXIS_NUMS)) integer() else values$TAXIS_NUMS,
             TAXIS_FLOATS = if (is.null(values$TAXIS_FLOATS)) numeric() else values$TAXIS_FLOATS,
             TAXIS_OFFSETS = if (is.null(values$TAXIS_OFFSETS)) numeric() else values$TAXIS_OFFSETS,
             IDCODE_STRING = if (is.null(values$IDCODE_STRING)) character() else values$IDCODE_STRING,
             IDCODE_DATE = if (is.null(values$IDCODE_DATE)) character() else values$IDCODE_DATE,
             BYTEORDER_STRING = if (is.null(values$BYTEORDER_STRING)) character() else values$BYTEORDER_STRING,
             BRICK_STATS = if (is.null(values$BRICK_STATS)) numeric() else values$BRICK_STATS,
             BRICK_TYPES = if (is.null(values$BRICK_TYPES)) integer() else values$BRICK_TYPES,
             BRICK_FLOAT_FACS = if (is.null(values$BRICK_FLOAT_FACS)) numeric() else values$BRICK_FLOAT_FACS,
             BRICK_LABS = if (is.null(values$BRICK_LABS)) character() else values$BRICK_LABS,
             BRICK_STATAUX = if (is.null(values$BRICK_STATAUX)) numeric() else values$BRICK_STATAUX,
             STAT_AUX = if (is.null(values$STAT_AUX)) numeric() else values$STAT_AUX,
             HISTORY_NOTE = if (is.null(values$HISTORY_NOTE)) character() else values$HISTORY_NOTE,
             NOTES_COUNT = if (is.null(values$NOTES_COUNT)) integer() else values$NOTES_COUNT,
             NOTE_NUMBER = if (is.null(values$NOTE_NUMBER)) character() else values$NOTE_NUMBER,
             TAGALIGN_MATVEC = if (is.null(values$TAGALIGN_MATVEC)) numeric() else values$TAGALIGN_MATVEC,
             VOLREG_MATVEC = if (is.null(values$VOLREG_MATVEC)) as.array(numeric()) else values$VOLREG_MATVEC,
             VOLREG_ROTCOM = if (is.null(values$VOLREG_ROTCOM)) character() else values$VOLREG_ROTCOM,
             VOLREG_CENTER_OLD = if (is.null(values$VOLREG_CENTER_OLD)) numeric() else values$VOLREG_CENTER_OLD,
             VOLREG_CENTER_BASE = if (is.null(values$VOLREG_CENTER_BASE)) numeric() else values$VOLREG_CENTER_BASE,
             VOLREG_ROTPARENT_IDCODE = if (is.null(values$VOLREG_ROTPARENT_IDCODE)) character() else values$VOLREG_ROTPARENT_IDCODE,
             VOLREG_ROTPARENT_NAME = if (is.null(values$VOLREG_ROTPARENT_NAME)) character() else values$VOLREG_ROTPARENT_NAME,
             VOLREG_GRIDPARENT_IDCODE = if (is.null(values$VOLREG_GRIDPARENT_IDCODE)) character() else values$VOLREG_GRIDPARENT_IDCODE,
             VOLREG_GRIDPARENT_NAME = if (is.null(values$VOLREG_GRIDPARENT_NAME)) character() else values$VOLREG_GRIDPARENT_NAME,
             VOLREG_INPUT_IDCODE = if (is.null(values$VOLREG_INPUT_IDCODE)) character() else values$VOLREG_INPUT_IDCODE,
             VOLREG_INPUT_NAME = if (is.null(values$VOLREG_INPUT_NAME)) character() else values$VOLREG_INPUT_NAME,
             VOLREG_BASE_IDCODE = if (is.null(values$VOLREG_BASE_IDCODE)) character() else values$VOLREG_BASE_IDCODE,
             VOLREG_BASE_NAME = if (is.null(values$VOLREG_BASE_NAME)) character() else values$VOLREG_BASE_NAME,
             VOLREG_ROTCOM_NUM = if (is.null(values$VOLREG_ROTCOM_NUM)) integer() else values$VOLREG_ROTCOM_NUM,
             IDCODE_ANAT_PARENT = if (is.null(values$IDCODE_ANAT_PARENT)) character() else values$IDCODE_ANAT_PARENT,
             TO3D_ZPAD = if (is.null(values$TO3D_ZPAD)) integer() else values$TO3D_ZPAD,
             IDCODE_WARP_PARENT = if (is.null(values$IDCODE_WARP_PARENT)) character() else values$IDCODE_WARP_PARENT,
             WARP_TYPE = if (is.null(values$WARP_TYPE)) integer() else values$WARP_TYPE,
             WARP_DATA = if (is.null(values$WARP_DATA)) numeric() else values$WARP_DATA,
             MARKS_XYZ = if (is.null(values$MARKS_XYZ)) numeric() else values$MARKS_XYZ,
             MARKS_LAB = if (is.null(values$MARKS_LAB)) character() else values$MARKS_LAB,
             MARKS_HELP = if (is.null(values$MARKS_HELP)) character() else values$MARKS_HELP,
             MARKS_FLAGS = if (is.null(values$MARKS_FLAGS)) integer() else values$MARKS_FLAGS,
             TAGSET_NUM = if (is.null(values$TAGSET_NUM)) integer() else values$TAGSET_NUM,
             TAGSET_FLOATS = if (is.null(values$TAGSET_FLOATS)) numeric() else values$TAGSET_FLOATS,
             TAGSET_LABELS = if (is.null(values$TAGSET_LABELS)) character() else values$TAGSET_LABELS,
             LABEL_1 = if (is.null(values$LABEL_1)) character() else values$LABEL_1,
             LABEL_2 = if (is.null(values$LABEL_2)) character() else values$LABEL_2,
             DATASET_NAME = if (is.null(values$DATASET_NAME)) character() else values$DATASET_NAME,
             DATASET_KEYWORDS = if (is.null(values$DATASET_KEYWORDS)) character() else values$DATASET_KEYWORDS,
             BRICK_KEYWORDS = if (is.null(values$BRICK_KEYWORDS)) character() else values$BRICK_KEYWORDS
             )
  ## Warnings?
  options(warn=oldwarn)
  return(nim)
}

############################################################################
## for writeS4.R
############################################################################
  
.AFNIheaderpart <- function(name, value, conhead) {
  a <- "\n"
  type <- switch(typeof(value),
                 "integer" = "integer",
                 "character" = "string",
                 "double" = "float")
  a <- paste(a, "type = ", type, "-attribute\n", sep="")
  a <- paste(a, "name = ", name, "\n", sep="")

  if (regexpr("string", type) == 1) {
    value <- paste("'", value, "~", sep="")
    a <- paste(a, "count = ", nchar(value) - 1, "\n", sep ="")
    a <- paste(a, value, "\n", sep="")
  } else {
    a <- paste(a, "count = ", length(value), "\n", sep ="")
    j <- 0
    while (j < length(value)) {
      left <- length(value) - j
      if (left > 4) {
        left <- 5
      }
      a <- paste(a, paste(value[(j+1):(j+left)], collapse="  "), "\n",
                 sep="  ")
      j <- j + 5
    }
  }
  writeChar(a, conhead, eos=NULL)
}

#' @title writeAFNI
#'   
#' @description This function saves a afni-class object to HEAD/BRIK pair in
#' AFNI format.
#'   
#' @details The \code{writeAFNI} function utilizes the internal \code{writeBin}
#' and \code{writeLines} command to write information to header/binary file
#' pair.
#'   
#' Current acceptable data types include 
#' \describe{ 
#' \item{INT16"}{DT SIGNED SHORT (16 bits per voxel)} 
#' \item{FLOAT32"}{DT FLOAT (32 bits per voxel)} 
#' \item{"COMPLEX128"}{DT COMPLEX (128 bits per voxel)} 
#' }
#'   
#' @name writeAFNI-methods
#' @aliases writeAFNI writeAFNI-methods writeAFNI,afni-method writeAFNI,ANY-method
#' @docType methods
#' @param nim is an object of class \code{afni}.
#' @param ... Additional variables defined by the method.
#' @param fname is the path and file name to save the AFNI file (.HEAD/BRIK)
#' \bold{without} the suffix.
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @param warn is a number to regulate the display of warnings (default = -1).
#' See \code{\link{options}} for more details.
#' @return Nothing.
#' @section Methods: 
#' \describe{ 
#' \item{nim = "afni"}{Write AFNI volume to disk.} 
#' \item{nim = "ANY"}{Not implemented.} 
#' }
#' @author Karsten Tabelow \email{karsten.tabelow@@wias-berlin.de}
#' @seealso \code{\link{writeANALYZE}}, \code{\link{writeNIfTI}}
#' @references AFNI\cr
#' \url{http://afni.nimh.nih.gov/pub/dist/src/README.attributes}
#' @keywords file methods
#' @examples
#' 
#' ## Taken from the AFNI Matlab Library
#' ## http://afni.nimh.nih.gov/pub/dist/data/afni_matlab_data.tgz
#' afni.path <- system.file("afni", package="oro.nifti")
#' orig <- readAFNI(file.path(afni.path, "ARzs_CW_avvr.DEL+orig"))
#' writeAFNI(orig, "test-afni-image", verbose=TRUE)
#' 
#' data <- readAFNI("test-afni-image", verbose=TRUE)
#' image(orig, zlim=c(0.5,256), oma=rep(2,4), bg="white")
#' image(data, zlim=c(0.5,256), oma=rep(2,4), bg="white")
#' abs.err <- abs(data - orig)
#' image(as(abs.err, "nifti"), zlim=range(0,1), oma=rep(2,4),
#'       bg="white")
#' @export
#' @rdname write_afni-methods
#' @docType methods
setGeneric("writeAFNI", function(nim,  ...) standardGeneric("writeAFNI"))
#' @export
#' @rdname write_afni-methods
#' @aliases writeAFNI,afni-method
setMethod("writeAFNI", "afni", function(nim, fname, verbose=FALSE, warn=-1) {
  ## Warnings?
  oldwarn <- getOption("warn")
  options(warn=warn)
  ## checking dataset!!!
  validAFNI <- getValidity(getClassDef("afni"))
  validAFNI(nim)
  ddim <- dim(nim@.Data)
  if (length(ddim) == 3) {
    dim(nim@.Data) <- ddim <- c(ddim, 1)
  }
  ## open header file for writing
  conhead <- file(paste(fname, ".HEAD", sep=""), "w")
  ## write mandatory attributes
  .AFNIheaderpart("DATASET_RANK", nim@DATASET_RANK, conhead)
  .AFNIheaderpart("DATASET_DIMENSIONS", nim@DATASET_DIMENSIONS, conhead)
  .AFNIheaderpart("TYPESTRING", nim@TYPESTRING, conhead)
  .AFNIheaderpart("SCENE_DATA", nim@SCENE_DATA, conhead)
  .AFNIheaderpart("ORIENT_SPECIFIC", nim@ORIENT_SPECIFIC, conhead)
  .AFNIheaderpart("ORIGIN", nim@ORIGIN, conhead)
  .AFNIheaderpart("DELTA", nim@DELTA, conhead)
  ## write mandatory attributes for time series
  if (length(nim@TAXIS_NUMS) > 0)
    .AFNIheaderpart("TAXIS_NUMS", nim@TAXIS_NUMS, conhead)
  if (length(nim@TAXIS_FLOATS) > 0)
    .AFNIheaderpart("TAXIS_FLOATS", nim@TAXIS_FLOATS, conhead)
  if (length(nim@TAXIS_OFFSETS) > 0)
    .AFNIheaderpart("TAXIS_OFFSETS", nim@TAXIS_OFFSETS, conhead)
  ## almost mandatory attributes
  if (length(nim@IDCODE_STRING) > 0) {
    .AFNIheaderpart("IDCODE_STRING", nim@IDCODE_STRING, conhead)
  } else {
    .AFNIheaderpart("IDCODE_DATE", "no-id", conhead)
  }
  if (length(nim@IDCODE_DATE) > 0) {
    .AFNIheaderpart("IDCODE_DATE", nim@IDCODE_DATE, conhead)
  } else {
    .AFNIheaderpart("IDCODE_DATE", date(), conhead)
  }
  if (length(nim@BYTEORDER_STRING) > 0) {
    endian <- switch(nim@BYTEORDER_STRING,
                     "MSB_FIRST" = "big",
                     "LSB_FIRST" = "little")
    .AFNIheaderpart("BYTEORDER_STRING", nim@BYTEORDER_STRING, conhead)
  } else {
    endian <- .Platform$endian
    .AFNIheaderpart("BYTEORDER_STRING", switch(endian,
                                               "little" = "LSB_FIRST",
                                               "big" = "MSB_FIRST"), conhead)
  }
  if (length(nim@BRICK_STATS) > 0) {
    .AFNIheaderpart("BRICK_STATS", nim@BRICK_STATS, conhead)
  } else {
    .AFNIheaderpart("BRICK_STATS", apply(nim@.Data, 4, range), conhead)
  }
  .AFNIheaderpart("BRICK_TYPES", nim@BRICK_TYPES, conhead)
  BRICK_FLOAT_FACS <- nim@BRICK_FLOAT_FACS
  if ((nim@BRICK_TYPES[1] == 1) && (max(abs(nim@BRICK_STATS)) > 32767)) { 
    for (k in 1:ddim[4]) {
      BRICK_FLOAT_FACS[k] <- max(abs(nim@BRICK_STATS[2*k-1]),
                                 abs(nim@BRICK_STATS[2*k])) / 32767
      nim@.Data[,,,k] <- nim@.Data[,,,k] / nim@BRICK_FLOAT_FACS[k]
    }
  }
  if (length(BRICK_FLOAT_FACS) > 0)
    .AFNIheaderpart("BRICK_FLOAT_FACS", BRICK_FLOAT_FACS, conhead)
  if (length(nim@BRICK_LABS) > 0)
    .AFNIheaderpart("BRICK_LABS", nim@BRICK_LABS, conhead)
  if (length(nim@BRICK_STATAUX) > 0)
    .AFNIheaderpart("BRICK_STATAUX", nim@BRICK_STATAUX, conhead)
  if (length(nim@STAT_AUX) > 0)
    .AFNIheaderpart("STAT_AUX", nim@STAT_AUX, conhead)
  ## note attributes
  if (length(nim@HISTORY_NOTE) > 0)
    .AFNIheaderpart("HISTORY_NOTE", nim@HISTORY_NOTE, conhead)
  if (length(nim@NOTES_COUNT) > 0)
    .AFNIheaderpart("NOTES_COUNT", nim@NOTES_COUNT, conhead)
  if (length(nim@NOTE_NUMBER) > 0) {
    for (i in 1:length(nim@NOTE_NUMBER)) {
      .AFNIheaderpart(paste(c("NOTE_NUMBER_",
                              rep("0", 3-nchar(as.character(i))),
                              as.character(i)), collapse=""),
                      nim@NOTE_NUMBER[i], conhead)
    }
  }
  ## Registration Attributes
  if (length(nim@TAGALIGN_MATVEC) > 0)
    .AFNIheaderpart("TAGALIGN_MATVEC", nim@TAGALIGN_MATVEC, conhead)
  if (length(nim@VOLREG_MATVEC) > 0) {
    for (i in 1:dim(nim@VOLREG_MATVEC)[2]) {
      .AFNIheaderpart(paste(c("VOLREG_MATVEC_",
                              rep("0", 6-nchar(as.character(i))),
                              as.character(i)), collapse=""),
                      as.numeric(nim@VOLREG_MATVEC[,i]), conhead)
    }
  }
  if (length(nim@VOLREG_ROTCOM) > 0) {
    for (i in 1:length(nim@VOLREG_ROTCOM)) {
      .AFNIheaderpart(paste(c("VOLREG_ROTCOM_",
                              rep("0", 6-nchar(as.character(i))),
                              as.character(i)), collapse=""),
                      nim@VOLREG_ROTCOM[i], conhead)
    }
  }
  if (length(nim@VOLREG_CENTER_OLD) > 0)
    .AFNIheaderpart("VOLREG_CENTER_OLD", nim@VOLREG_CENTER_OLD, conhead)
  if (length(nim@VOLREG_CENTER_BASE) > 0)
    .AFNIheaderpart("VOLREG_CENTER_BASE", nim@VOLREG_CENTER_BASE, conhead)
  if (length(nim@VOLREG_ROTPARENT_IDCODE) > 0)
    .AFNIheaderpart("VOLREG_ROTPARENT_IDCODE", nim@VOLREG_ROTPARENT_IDCODE, conhead)
  if (length(nim@VOLREG_ROTPARENT_NAME) > 0)
    .AFNIheaderpart("VOLREG_ROTPARENT_NAME", nim@VOLREG_ROTPARENT_NAME, conhead)
  if (length(nim@VOLREG_GRIDPARENT_IDCODE) > 0)
    .AFNIheaderpart("VOLREG_GRIDPARENT_IDCODE", nim@VOLREG_GRIDPARENT_IDCODE, conhead)
  if (length(nim@VOLREG_GRIDPARENT_NAME) > 0)
    .AFNIheaderpart("VOLREG_GRIDPARENT_NAME", nim@VOLREG_GRIDPARENT_NAME, conhead)
  if (length(nim@VOLREG_INPUT_IDCODE) > 0)
    .AFNIheaderpart("VOLREG_INPUT_IDCODE", nim@VOLREG_INPUT_IDCODE, conhead)
  if (length(nim@VOLREG_INPUT_NAME) > 0)
    .AFNIheaderpart("VOLREG_INPUT_NAME", nim@VOLREG_INPUT_NAME, conhead)
  if (length(nim@VOLREG_BASE_IDCODE) > 0)
    .AFNIheaderpart("VOLREG_BASE_IDCODE", nim@VOLREG_BASE_IDCODE, conhead)
  if (length(nim@VOLREG_BASE_NAME) > 0)
    .AFNIheaderpart("VOLREG_BASE_NAME", nim@VOLREG_BASE_NAME, conhead)
  if (length(nim@VOLREG_ROTCOM_NUM) > 0)
    .AFNIheaderpart("VOLREG_ROTCOM_NUM", nim@VOLREG_ROTCOM_NUM, conhead)

  ## Miscellaneous Attributes
  if (length(nim@IDCODE_ANAT_PARENT) > 0)
    .AFNIheaderpart("IDCODE_ANAT_PARENT", nim@IDCODE_ANAT_PARENT, conhead)
  if (length(nim@TO3D_ZPAD) > 0)
    .AFNIheaderpart("TO3D_ZPAD", nim@TO3D_ZPAD, conhead)

  ## Warping Attributes
  if (length(nim@IDCODE_WARP_PARENT) > 0)
    .AFNIheaderpart("IDCODE_WARP_PARENT", nim@IDCODE_WARP_PARENT, conhead)
  if (length(nim@WARP_TYPE) > 0)
    .AFNIheaderpart("WARP_TYPE", nim@WARP_TYPE, conhead)
  if (length(nim@WARP_DATA) > 0)
    .AFNIheaderpart("WARP_DATA", nim@WARP_DATA, conhead)

  ## Talairach Markers Attributes
  if (length(nim@MARKS_XYZ) > 0)
    .AFNIheaderpart("MARKS_XYZ", nim@MARKS_XYZ, conhead)
  if (length(nim@MARKS_LAB) > 0)
    .AFNIheaderpart("MARKS_LAB", nim@MARKS_LAB, conhead)
  if (length(nim@MARKS_HELP) > 0)
    .AFNIheaderpart("MARKS_HELP", nim@MARKS_HELP, conhead)
  if (length(nim@MARKS_FLAGS) > 0)
    .AFNIheaderpart("MARKS_FLAGS", nim@MARKS_FLAGS, conhead)

  ## Attributes for User-Defined Tags
  if (length(nim@TAGSET_NUM) > 0)
    .AFNIheaderpart("TAGSET_NUM", nim@TAGSET_NUM, conhead)
  if (length(nim@TAGSET_FLOATS) > 0)
    .AFNIheaderpart("TAGSET_FLOATS", nim@TAGSET_FLOATS, conhead)
  if (length(nim@TAGSET_LABELS) > 0)
    .AFNIheaderpart("TAGSET_LABELS", nim@TAGSET_LABELS, conhead)

  ## Nearly Useless Attributes
  if (length(nim@LABEL_1) > 0)
    .AFNIheaderpart("LABEL_1", nim@LABEL_1, conhead)
  if (length(nim@LABEL_2) > 0)
    .AFNIheaderpart("LABEL_2", nim@LABEL_2, conhead)
  if (length(nim@DATASET_NAME) > 0)
    .AFNIheaderpart("DATASET_NAME", nim@DATASET_NAME, conhead)
  if (length(nim@DATASET_KEYWORDS) > 0)
    .AFNIheaderpart("DATASET_KEYWORDS", nim@DATASET_KEYWORDS, conhead)
  if (length(nim@BRICK_KEYWORDS) > 0)
    .AFNIheaderpart("BRICK_KEYWORDS", nim@BRICK_KEYWORDS, conhead)
  
  ## close header file
  close(conhead)
  ## write data
  if (!(nim@BRICK_TYPES[1] %in% c(1,3,5))) {
    stop("Sorry, cannot write this BRICK_TYPES.", call.=FALSE)
  }
  conbrik <- file(paste(fname, ".BRIK", sep=""), "wb")
  switch(as.character(nim@BRICK_TYPES[1]),
         "1" = writeBin(as.integer(nim@.Data), conbrik, size=2, endian=endian),
         "3" = writeBin(as.numeric(nim@.Data), conbrik, size=4, endian=endian),
         "5" = writeBin(as.complex(nim@.Data), conbrik, size=16, endian=endian))
  close(conbrik)
  ## Warnings?
  options(warn=oldwarn)
  invisible()
})

##############################################################################
## for convert_afni.R
##############################################################################

#' @title Convert AFNI data codes
#' 
#' @description Codes that appear in the AFNI header are mapped to meaningful character
#' strings.
#' 
#' @details \code{switch} statements are used to map a numeric code to the appropriate
#' string.
#' 
#' @name convert.scene
#' @docType data
#' @param scene.data defines data type.
#' @param typestring defines whether func or anat data.
#' @return A character string.
#' @author Karsten Tabelow \email{karsten.tabelow@@wias-berlin.de}
#' @seealso \code{\link{convert.datatype.anlz}},
#' \code{\link{convert.orient.anlz}}
#' @references AFNI\cr
#' \url{http://afni.nimh.nih.gov/pub/dist/src/README.attributes}
#' @keywords misc
#' @examples
#' 
#' ## 4 = CT for anatomic data
#' convert.scene(4, "3DIM_HEAD_ANAT")
#' @export 
#' @rdname convert_afni
convert.scene <- function (scene.data, typestring) {
  if (grepl("ANAT", typestring)) {
    switch(as.character(scene.data),
           "0" = "SPGR",
           "1" = "FSE",
           "2" = "EPI",
           "3" = "MRAN",
           "4" = "CT",
           "5" = "SPECT",
           "6" = "PET",
           "7" = "MRA",
           "8" = "BMAP",
           "9" = "DIFF",
           "10" = "OMRI",
           "11" = "BUCK")
  } else if (grepl("FUNC", typestring)) {
    switch(as.character(scene.data),
           "0" = "FIM",
           "1" = "THR",
           "2" = "COR",
           "3" = "TT",
           "4" = "FT",
           "5" = "ZT",
           "6" = "CT",
           "7" = "BT",
           "8" = "BN",
           "9" = "GT",
           "10" = "PT",
           "11" = "BUCK")  
  } else {
    ""
  }
}


