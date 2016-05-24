##
##
## Copyright (c) 2009-2014 Brandon Whitcher and Volker Schmid
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
## $Id: readS4.R 332 2010-01-29 16:54:07Z bjw34032 $
##

## Sub-routines
.readCharWithEmbeddedNuls <- function(fid, n, to="UTF-8") {
  txt <- readBin(fid, "raw", n)
  iconv(rawToChar(txt[txt != as.raw(0)]), to=to)
}

##
## readNIfTI() is a convient interface for the user
##
#'
#' @title readNIfTI
#' 
#' @description These functions read in the header information and multidimensional array
#' from a binary file in NIfTI-1 format into a \code{\linkS4class{nifti}}-class
#' object.
#' 
#' @details The \code{readNIfTI} function utilizes internal methods \code{readBin} and
#' \code{readChar} to efficiently extract information from the binary file(s).
#' 
#' Current acceptable data types include \describe{ \item{list("UINT8")}{BINARY
#' (1 bit per voxel)} \item{list("INT16")}{SIGNED SHORT (16 bits per voxel)}
#' \item{list("INT32")}{SINGED INT (32 bits per voxel)}
#' \item{list("FLOAT32")}{FLOAT (32 bits per voxel)}
#' \item{list("DOUBLE64")}{DOUBLE (64 bits per voxel)}
#' \item{list("UINT16")}{UNSIGNED SHORT (16 bits per voxel)}
#' \item{list("UINT32")}{UNSIGNED INT (32 bits per voxel)} }
#' 
#' @param fname is the file name of the NIfTI file(s).
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @param warn is a number to regulate the display of warnings (default = -1).
#' See \code{options} for more details.
#' @param reorient is a logical variable (default = \code{TRUE}) that enforces
#' Qform/Sform transformations.
#' @param call keeps track of the current function call for use in the NIfTI
#' extension.
#' @return An object of class \code{nifti}.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com},\cr Volker Schmid
#' \email{volkerschmid@@users.sourceforge.net},\cr Andrew Thornton
#' \email{zeripath@@users.sourceforge.net}
#' @seealso \code{\link{readAFNI}}, \code{\link{readANALYZE}}
#' @references NIfTI-1\cr \url{http://nifti.nimh.nih.gov/}
#' @keywords file
#' @examples
#' 
#' \dontrun{
#' url <- "http://nifti.nimh.nih.gov/nifti-1/data/filtered_func_data.nii.gz"
#' urlfile <- file.path(system.file("nifti", package="oro.nifti"),
#'                      "filtered_func_data")
#' download.file(url, urlfile, quiet=TRUE)
#' }
#' ## The NIfTI file provided here contains the first 18 volumes (10%)
#' ## of the original data set
#' urlfile <- file.path(system.file("nifti", package="oro.nifti"),
#'                      "filtered_func_data")
#' (ffd <- readNIfTI(urlfile))
#' image(ffd, oma=rep(2,4))
#' orthographic(ffd, oma=rep(2,4))
#' \dontrun{
#' ## 27 scans of Colin Holmes (MNI) brain co-registered and averaged
#' ## NIfTI two-file format
#' URL <- "http://imaging.mrc-cbu.cam.ac.uk/downloads/Colin/colin_1mm.tgz"
#' urlfile <- file.path(tempdir(), "colin_1mm.tgz")
#' download.file(URL, dest=urlfile, quiet=TRUE)
#' untar(urlfile, exdir=tempdir())
#' colin <- readNIfTI(file.path(tempdir(), "colin_1mm"))
#' image(colin, oma=rep(2,4))
#' orthographic(colin, oma=rep(2,4))
#' }
#' @rdname read_nifti
#' @export
#' @name readNIfTI
readNIfTI <- function(fname, verbose=FALSE, warn=-1, reorient=TRUE,
                      call=NULL) {
  if (is.null(call)) {
    call <- match.call()
  }
  ## Warnings?
  oldwarn <- getOption("warn")
  options(warn=warn)
  if (verbose) {
    cat(paste("  fname =", fname), fill=TRUE)
  }
  ## Separate path and file name for more robust extension stripping
  pathVector <- unlist(strsplit(fname, "/"))
  file.name <- pathVector[length(pathVector)]
  path <- paste(pathVector[-length(pathVector)], collapse="/")
  if (length(pathVector) > 1) {
      fname <- paste(path, file.name, sep="/")
  } else {
      fname <- file.name
  }
  ## Strip any extensions
  fname <- sub("\\.gz$", "", fname)
  fname <- sub("\\.nii$", "", fname)
  fname <- sub("\\.hdr$", "", fname)
  fname <- sub("\\.img$", "", fname)
  ## Add all possible file extensions 
  nii <- paste(fname, "nii", sep=".")
  niigz <- paste(fname, "nii.gz", sep=".")
  hdr <- paste(fname, "hdr", sep=".")
  hdrgz <- paste(fname, "hdr.gz", sep=".")
  img <- paste(fname, "img", sep=".")
  imggz <- paste(fname, "img.gz", sep=".")
  ## Check all possible file extensions
  if (file.exists(niigz)) {
    ## If compressed file exists, then upload!
    if (verbose) {
      cat(paste("  files =", niigz), fill=TRUE)
    }
    nim <- .read.nifti.content(fname, onefile=TRUE, gzipped=TRUE,
                               verbose=verbose, warn=warn, reorient=reorient,
                               call=call)
  } else {
    if (file.exists(nii)) {
      ## If uncompressed file exists, then upload!
      if (verbose) {
        cat(paste("  files =", nii), fill=TRUE)
      }
      nim <- .read.nifti.content(fname, onefile=TRUE, gzipped=FALSE,
                                 verbose=verbose, warn=warn, reorient=reorient,
                                 call=call)
    } else {
      if (file.exists(hdrgz) && file.exists(imggz)) {
        ## If compressed files exist, then upload!
        if (verbose) {
          cat(paste("  files =", hdrgz, "and", imggz), fill=TRUE)
        }
        nim <- .read.nifti.content(fname, onefile=FALSE, gzipped=TRUE,
                                   verbose=verbose, warn=warn,
                                   reorient=reorient, call=call)
      } else {
        ## If uncompressed files exist, then upload!
        if (file.exists(hdr) && file.exists(img)) {
          if (verbose) {
            cat(paste("  files =", hdr, "and", img), fill=TRUE)
          }
        nim <- .read.nifti.content(fname, onefile=FALSE, gzipped=FALSE,
                                   verbose=verbose, warn=warn,
                                   reorient=reorient, call=call)
        } else {
          stop("File(s) not found!")
        }
      }
    }
  }
  ### Reset cal_max and cal_min - in case these do not work correctly
  nim = calibrateImage(nim, infok = TRUE)
  options(warn=oldwarn)
  return(nim)
}

############################################################################
############################################################################
############################################################################

.read.nifti.content <- function(fname, onefile=TRUE, gzipped=TRUE,
                                verbose=FALSE, warn=-1, reorient=FALSE,
                                call=NULL) {
  ## Open appropriate file
  if (gzipped) {
    suffix <- ifelse(onefile, "nii.gz", "hdr.gz")
    fname <- paste(fname, suffix, sep=".")
    fid <- gzfile(fname, "rb")
    if (verbose) {
      cat("  nii   =", fname, fill=TRUE)
    }
  } else {
    suffix <- ifelse(onefile, "nii", "hdr")
    fname <- paste(fname, suffix, sep=".")
    fid <- file(fname, "rb")
    if (verbose) {
      cat("  hdr   =", fname, fill=TRUE)
    }
  }
  ## Warnings?
  oldwarn <- getOption("warn")
  options(warn=warn)
  ## Test for endian properties
  endian <- .Platform$endian
  sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
  if (sizeof.hdr != 348) {
    close(fid)
    endian <- "swap"
    if (gzipped) {
      fid <- gzfile(fname, "rb")
    } else {
      fid <- file(fname, "rb")
    }
    sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
    if (verbose) {
      cat("  ENDIAN = swap", fill=TRUE)
    }
  }
  ## Construct S4 object
  nim <- nifti()
  nim@"sizeof_hdr" <- sizeof.hdr
  nim@"data_type" <- .readCharWithEmbeddedNuls(fid, n=10)
  nim@"db_name" <- .readCharWithEmbeddedNuls(fid, n=18)
  nim@"extents" <- readBin(fid, integer(), size=4, endian=endian)
  nim@"session_error" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"regular" <- .readCharWithEmbeddedNuls(fid, n=1)
  nim@"dim_info" <- .readCharWithEmbeddedNuls(fid, n=1)
  nim@"dim_" <- readBin(fid, integer(), 8, size=2, endian=endian)
  nim@"intent_p1" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"intent_p2" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"intent_p3" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"intent_code" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"datatype" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"bitpix" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"slice_start" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"pixdim" <- readBin(fid, numeric(), 8, size=4, endian=endian)
  nim@"vox_offset" <- readBin(fid, numeric(), size=4, endian=endian)
  if (verbose) {
    cat("  vox_offset =", nim@"vox_offset", fill=TRUE)
  }
  nim@"scl_slope" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"scl_inter" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"slice_end" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"slice_code" <- readBin(fid, integer(), size=1, signed=FALSE,
                              endian=endian)
  nim@"xyzt_units" <- readBin(fid, integer(), size=1, signed=FALSE,
                              endian=endian)
  nim@"cal_max" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"cal_min" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"slice_duration" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"toffset" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"glmax" <- readBin(fid, integer(), size=4, endian=endian)
  nim@"glmin" <- readBin(fid, integer(), size=4, endian=endian)
  nim@"descrip" <- .readCharWithEmbeddedNuls(fid, n=80)
  nim@"aux_file" <- .readCharWithEmbeddedNuls(fid, n=24)
  nim@"qform_code" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"sform_code" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"quatern_b" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"quatern_c" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"quatern_d" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"qoffset_x" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"qoffset_y" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"qoffset_z" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"srow_x" <- readBin(fid, numeric(), 4, size=4, endian=endian)
  nim@"srow_y" <- readBin(fid, numeric(), 4, size=4, endian=endian)
  nim@"srow_z" <- readBin(fid, numeric(), 4, size=4, endian=endian)
  nim@"intent_name" <- .readCharWithEmbeddedNuls(fid, n=16)
  nim@"magic" <- .readCharWithEmbeddedNuls(fid, n=4)
  ## To flag such a struct as being conformant to the NIFTI-1 spec,
  ## the last 4 bytes of the header must be either the C String "ni1"
  ## or "n+1"; in hexadecimal, the 4 bytes [6E 69 31 00] or [6E 2B 31
  ## 00] (in any future version of this format, the 1 will be upgraded
  ## to 2, etc.).  Normally, such a "magic number" or flag goes at the
  ## start of the file, but trying to avoid clobbering widely-used
  ## ANALYZE 7.5 fields led to putting this marker last.  However,
  ## recall that "the last shall be first" (Matthew 20:16).
  if (onefile) {  
    if (nim@"magic" != "n+1") {
      stop("This is not a one-file NIfTI format")
    }   
    nim@"extender" <- readBin(fid, integer(), 4, size=1, signed=FALSE,
                              endian=endian)
    ## If extension[0] is nonzero, it indicates that extended header
    ## information is present in the bytes following the extension
    ## array.  In a .nii file, this extended header data is before the
    ## image data (and vox_offset must be set correctly to allow for
    ## this).  In a .hdr file, this extended data follows extension and
    ## proceeds (potentially) to the end of the file.
    ##
    if (nim@"extender"[1] > 0 || nim@"vox_offset" > 352) {
      if (verbose) {
        cat("  niftiExtension detected!", fill=TRUE)
      }
      if (!is(nim, "niftiExtension")) {
        nim <- as(nim, "niftiExtension")
      }
      while (seek(fid) < nim@"vox_offset") {
        if (verbose) {
          cat("  seek(fid) =", seek(fid), fill=TRUE)
        }
        nimextsec <- new("niftiExtensionSection")
        nimextsec@esize <- readBin(fid, integer(), size=4, endian=endian)
        nimextsec@ecode <- readBin(fid, integer(), size=4, endian=endian)
        nimextsec@edata <- .readCharWithEmbeddedNuls(fid, n=nimextsec@esize-8)
        nim@extensions <- append(nim@extensions, nimextsec)
      }
      if (seek(fid) > nim@"vox_offset") {
        stop("-- extension size (esize) has overshot voxel offset --")
      }
    }
  }

  if (verbose) {
    cat("  seek(fid) =", seek(fid), fill=TRUE)
  }
  dims <- 2:(1+nim@"dim_"[1])
  n <- prod(nim@"dim_"[dims])
  if (! onefile) {
    if (nim@"magic" != "ni1") {
      stop("This is not a two-file NIfTI format")
    }
    close(fid)
    fname <- sub("\\.hdr$", "\\.img", fname)
    if (gzipped) {
      fid <- gzfile(fname, "rb")
    } else {
      fid <- file(fname, "rb")
    }
    ## seek(fid, nim@"vox_offset") # not necessary for two-file format
  }
  data <-
    switch(as.character(nim@"datatype"),
           "2" = readBin(fid, integer(), n, nim@"bitpix"/8, signed=FALSE,
             endian=endian),
           "4" = readBin(fid, integer(), n, nim@"bitpix"/8, endian=endian),
           "8" = readBin(fid, integer(), n, nim@"bitpix"/8, endian=endian),
           "16" = readBin(fid, double(), n, nim@"bitpix"/8, endian=endian),
           "64" = readBin(fid, double(), n, nim@"bitpix"/8, endian=endian),
           "512" = readBin(fid, integer(), n, nim@"bitpix"/8, signed=FALSE,
             endian=endian),
           "768" = readBin(fid, integer(), n, nim@"bitpix"/8, signed=FALSE,
             endian=endian),
           stop(paste("Data type", nim@"datatype", "unsupported in", fname))
           )
  close(fid)

  ## WARNING to the user
  if (nim@"scl_slope" != 0) {
    warning(paste("scl_slope =", nim@"scl_slope", "and data must be rescaled."))
    data <- data * nim@"scl_slope" + nim@"scl_inter"
  }
  ##
  ## THE SLOW BIT FOLLOWS
  ##
  ## 3D IMAGE (VOLUME) ORIENTATION AND LOCATION IN SPACE:
  ## There are 3 different methods by which continuous coordinates can
  ## attached to voxels.  The discussion below emphasizes 3D volumes,
  ## and the continuous coordinates are referred to as (x,y,z).  The
  ## voxel index coordinates (i.e., the array indexes) are referred to
  ## as (i,j,k), with valid ranges:
  ##   i = 0 .. dim[1]-1
  ##   j = 0 .. dim[2]-1  (if dim[0] >= 2)
  ##   k = 0 .. dim[3]-1  (if dim[0] >= 3)
  ## The (x,y,z) coordinates refer to the CENTER of a voxel.  In
  ## methods 2 and 3, the (x,y,z) axes refer to a subject-based
  ## coordinate system, with
  ##   +x = Right  +y = Anterior  +z = Superior.
  ## This is a right-handed coordinate system.  However, the exact
  ## direction these axes point with respect to the subject depends on
  ## qform_code (Method 2) and sform_code (Method 3).
  ##
  if (reorient) {
    nim@.Data <- reorient(nim, data, verbose=verbose)
    nim@"reoriented" <- TRUE
  } else {
    nim@.Data <- array(data, nim@"dim_"[dims])
  }
  ## Warnings?
  options(warn=oldwarn)
  ## Check validity
  validNIfTI <- getValidity(getClassDef("nifti"))
  validNIfTI(nim)
  if (getOption("niftiAuditTrail")) {
    if (is.null(call)) {
      call <- match.call()
    }
    nim <- niftiExtensionToAuditTrail(nim, workingDirectory=getwd(),
                                      filename=fname, call=call)
  }
  return(nim)
}

############################################################################
## readANALYZE() is a convenient interface for the user
############################################################################
#' @title readANALYZE
#' 
#' @description These functions read in the header information and multi-dimensional array
#' from a binary file in Analyze 7.5 format.
#' 
#' @details The internal functions \code{readBin} and \code{rawToChar} are utilized in
#' order to efficiently extract information from a binary file.  The types of
#' data are limited to 1- and 2-byte integers, 4-byte floats and 8-byte
#' doubles.
#' 
#' @param fname Pathname of the Analyze pair of files .img and .hdr without the
#' suffix.
#' @param SPM is a logical variable (default = \code{FALSE}) that forces the
#' voxel data values to be rescaled using the funused1 ANALYZE header field.
#' This is an undocumented convention of ANALYZE files processed using the
#' Statistical Parametric Mapping (SPM) software.
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @param warn is a number to regulate the display of warnings (default = -1).
#' See \code{options} for more details.
#' @return An object of class \code{anlz} is produced.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com},\cr Volker Schmid
#' \email{volkerschmid@@users.sourceforge.net}
#' @seealso \code{\link{readNIfTI}}
#' @references ANALYZE 7.5\cr \url{https://rportal.mayo.edu/bir/ANALYZE75.pdf}
#' @keywords file
#' @examples
#' 
#' ## avg152T1
#' anlz.path <- system.file("anlz", package="oro.nifti")
#' mni152 <- readANALYZE(file.path(anlz.path, "avg152T1"))
#' image(mni152, oma=rep(2,4))
#' orthographic(mni152, oma=rep(2,4))
#' @name readANALYZE
#' @rdname read_anlz
#' @export
readANALYZE <- function(fname, SPM=FALSE, verbose=FALSE, warn=-1) {
  ## Warnings?
  oldwarn <- getOption("warn")
  options(warn=warn)
  if (verbose) {
    cat(paste("  fname =", fname), fill=TRUE)
  }
  ## Separate path and file name for more robust extension stripping
  pathVector <- unlist(strsplit(fname, "/"))
  file.name <- pathVector[length(pathVector)]
  path <- paste(pathVector[-length(pathVector)], collapse="/")
  if (length(pathVector) > 1) {
      fname <- paste(path, file.name, sep="/")
  } else {
      fname <- file.name
  }
  ## Strip any extensions
  fname <- sub("\\.gz$", "", fname)
  fname <- sub("\\.hdr$", "", fname)
  fname <- sub("\\.img$", "", fname)
  if (! (file.exists(paste(fname, "hdr", sep=".")) &&
         file.exists(paste(fname, "img", sep="."))) &&
      ! (file.exists(paste(fname, "hdr.gz", sep=".")) &&
         file.exists(paste(fname, "img.gz", sep=".")))) {
    stop("File(s) not found!")
  }
  ## If uncompressed files exist, then upload!
  if (file.exists(paste(fname, "hdr", sep=".")) &&
      file.exists(paste(fname, "img", sep="."))) {      
    if (verbose) {
      cat(paste("  files = ", fname, ".{hdr,img}", sep=""), fill=TRUE)
    }
    aim <- .read.analyze.content(fname, gzipped=FALSE, SPM=SPM,
                                 verbose=verbose, warn=warn)
    options(warn=oldwarn)
    return(aim)
  }
  ## If compressed files exist, then upload!
  if (file.exists(paste(fname, "hdr.gz", sep=".")) &&
      file.exists(paste(fname, "img.gz", sep="."))) {      
    if (verbose) {
      cat(paste("  files = ", fname, ".{hdr.gz,img.gz}", sep=""), fill=TRUE)
    }
    aim <- .read.analyze.content(fname, gzipped=TRUE, SPM=SPM,
                                 verbose=verbose, warn=warn)
    options(warn=oldwarn)
    return(aim)
  }
  invisible()
}

############################################################################
############################################################################
############################################################################

.read.analyze.content <- function(fname, gzipped=TRUE, SPM=FALSE,
                                  verbose=FALSE, warn=-1) {
  ## Open header file
  if (gzipped) {
    fname <- paste(fname, "hdr.gz", sep=".")
    fid <- gzfile(fname, "rb")
  } else {
    fname <- paste(fname, "hdr", sep=".")
    fid <- file(fname, "rb")
  }
  if (verbose) {
    cat("  hdr   =", fname, fill=TRUE)
  }
  ## Warnings?
  oldwarn <- getOption("warn")
  options(warn=warn)
  ## Test for endian properties
  endian <- .Platform$endian
  sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
  if (sizeof.hdr != 348) {
    close(fid)
    endian <- "swap"
    if (gzipped) {
      fid <- gzfile(fname, "rb")
    } else {
      fid <- file(fname, "rb")
    }
    sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
  }
  ## Construct S4 object
  aim <- new("anlz")
  aim@"sizeof_hdr" <- sizeof.hdr
  aim@"data_type" <- .readCharWithEmbeddedNuls(fid, n=10)
  aim@"db_name" <- .readCharWithEmbeddedNuls(fid, n=18)
  aim@"extents" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"session_error" <- readBin(fid, integer(), size=2, endian=endian)
  aim@"regular" <- .readCharWithEmbeddedNuls(fid, n=1)
  aim@"hkey_un0" <- .readCharWithEmbeddedNuls(fid, n=1)
  aim@"dim_" <- readBin(fid, integer(), 8, size=2, endian=endian)
  aim@"vox_units" <- .readCharWithEmbeddedNuls(fid, n=4)
  aim@"cal_units" <- .readCharWithEmbeddedNuls(fid, n=8)
  aim@"unused1" <- readBin(fid, integer(), size=2, endian=endian)
  aim@"datatype" <- readBin(fid, integer(), size=2, endian=endian)
  aim@"bitpix" <- readBin(fid, integer(), size=2, endian=endian)
  aim@"dim_un0" <- readBin(fid, integer(), size=2, endian=endian)
  aim@"pixdim" <- readBin(fid, numeric(), 8, size=4, endian=endian)
  aim@"vox_offset" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"funused1" <- readBin(fid, numeric(), size=4, endian=endian)
  ## SPM has used the ANALYZE 7.5 funused1 field as a scaling factor
  if (aim@"funused1" != 0) {
    warning(paste("funused1 =", aim@"funused1", "and data must be rescaled."))
  }
  ##
  aim@"funused2" <- readBin(fid, numeric(), size=4, endian=endian)
  ## Maybe I'm paranoid, but let's check funused2 in case it is an intercept
  if (aim@"funused2" != 0) {
    warning(paste("funused2 =", aim@"funused2", "and data must be translated."))
  }
  ##
  aim@"funused3" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"cal_max" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"cal_min" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"compressed" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"verified" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"glmax" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"glmin" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"descrip" <- .readCharWithEmbeddedNuls(fid, n=80)
  aim@"aux_file" <- .readCharWithEmbeddedNuls(fid, n=24)
  aim@"orient" <- .readCharWithEmbeddedNuls(fid, n=1)
  aim@"origin" <- readBin(fid, integer(), 5, size=2, endian=endian) # .readCharWithEmbeddedNuls(fid, n=10)
  aim@"generated" <- .readCharWithEmbeddedNuls(fid, n=10)
  aim@"scannum" <- .readCharWithEmbeddedNuls(fid, n=10)
  aim@"patient_id" <- .readCharWithEmbeddedNuls(fid, n=10)
  aim@"exp_date" <- .readCharWithEmbeddedNuls(fid, n=10)
  aim@"exp_time" <- .readCharWithEmbeddedNuls(fid, n=10)
  aim@"hist_un0" <- .readCharWithEmbeddedNuls(fid, n=3)
  aim@"views" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"vols_added" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"start_field" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"field_skip" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"omax" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"omin" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"smax" <- readBin(fid, integer(), size=4, endian=endian)
  magic <- readBin(fid, raw(), n=4, endian=endian)
  aim@"smin" <- readBin(magic, integer(), size=4, endian=endian)
  ## Test for "magic" field (should not exist)
  if (rawToChar(magic) == "ni1") { # now its actually NIfTI two-file format
    stop("This is in two-file NIfTI format, please use readNIfTI")
  }
  close(fid)
  ## Open image file
  if (gzipped) {
    fname <- sub("\\.hdr\\.gz$", "\\.img\\.gz", fname) # paste(fname, ".img.gz", sep=".")
    fid <- gzfile(fname, "rb")
  } else {
    fname <- sub("\\.hdr$", "\\.img", fname) # paste(fname, "img", sep=".")
    fid <- file(fname, "rb")
  }
  if (verbose) {
    cat("  img   =", fname, fill=TRUE)
  }
  n <- prod(aim@"dim_"[2:5])
  data <- switch(as.character(aim@"datatype"),
                 "1" = readBin(fid, integer(), n, aim@"bitpix"/8,
                   signed=FALSE, endian=endian),
                 "2" = readBin(fid, integer(), n, aim@"bitpix"/8,
                   signed=FALSE, endian=endian),
                 "4" = readBin(fid, integer(), n, aim@"bitpix"/8,
                   endian=endian),
                 "8" = readBin(fid, integer(), n, aim@"bitpix"/8,
                   endian=endian),
                 "16" = readBin(fid, numeric(), n, aim@"bitpix"/8,
                   endian=endian),
                 "64" = readBin(fid, double(), n, aim@"bitpix"/8,
                   endian=endian),
                 stop(paste("Data type ", aim@"datatype", " (",
                            convert.datatype.anlz(aim@"datatype"), 
                            ") unsupported in", fname, sep="")))
  close(fid)
  dims <- 2:(1+aim@"dim_"[1])
  if (SPM) {
    if (verbose) {
      cat("  SPM format has been specified and data re-scaled.", fill=TRUE)
    }
    aim@.Data <- array(aim@"funused1" * data, aim@"dim_"[dims])
  } else {
    aim@.Data <- array(data, aim@"dim_"[dims])
  }
  ## Warnings?
  options(warn=oldwarn)
  ### Reset cal_max and cal_min - in case these do not work correctly
  aim <- calibrateImage(aim, infok = TRUE)  
  ## Check validity
  validANALYZE <- getValidity(getClassDef("anlz"))
  validANALYZE(aim)
  return(aim)
}
