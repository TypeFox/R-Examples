##
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
## $Id: writeS4.R 332 2010-01-29 16:54:07Z bjw34032 $
##

setGeneric("writeNIfTI", function(nim,  ...) standardGeneric("writeNIfTI"))
#' @title writeNIfTI
#' 
#' @description This function saves a NIfTI-class object to a single binary file in NIfTI
#' format.
#' 
#' @details The \code{writeNIfTI} function utilizes the internal \code{writeBin} and
#' \code{writeChar} command to write information to a binary file.
#' 
#' Current acceptable data types include \describe{ \item{list("UINT8")}{DT
#' BINARY (1 bit per voxel)} \item{list("INT16")}{DT SIGNED SHORT (16 bits per
#' voxel)} \item{list("INT32")}{DT SINGED INT (32 bits per voxel)}
#' \item{list("FLOAT32")}{DT FLOAT (32 bits per voxel)}
#' \item{list("DOUBLE64")}{DT DOUBLE (64 bits per voxel)}
#' \item{list("UINT16")}{DT UNSIGNED SHORT (16 bits per voxel)} }
#' 
#' @name writeNIfTI-methods
#' @aliases writeNIfTI writeNIfTI-methods writeNIfTI,nifti-method writeNIfTI,anlz-method
#' writeNIfTI,array-method
#' @docType methods
#' @param nim is an object of class \code{nifti} or \code{anlz}.
#' @param filename is the path and file name to save the NIfTI file (.nii)
#' \bold{without} the suffix.
#' @param onefile is a logical value that allows the scanning of single-file
#' (.nii) or dual-file format (.hdr and .img) NIfTI files (default =
#' \code{TRUE}).
#' @param gzipped is a character string that enables exportation of compressed
#' (.gz) files (default = \code{TRUE}).
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @param warn is a number to regulate the display of warnings (default = -1).
#' See \code{\link{options}} for more details.
#' @return Nothing.
#' @section Methods: \describe{ \item{object = "anlz"}{Convert ANALYZE object
#' to class \code{nifti} and write the NIfTI volume to disk.} \item{object =
#' "array"}{Convert array to class \code{nifti} and write the NIfTI volume to
#' disk.} \item{object = "nifti"}{Write NIfTI volume to disk.} }
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com},\cr Volker Schmid
#' \email{volkerschmid@@users.sourceforge.net}
#' @seealso \code{\link{writeAFNI}}, \code{\link{writeANALYZE}}
#' @references NIfTI-1\cr \url{http://nifti.nimh.nih.gov/}
#' @import methods
#' @import utils
#' @keywords file methods
#' @examples
#' 
#' norm <- dnorm(seq(-5, 5, length=32), sd=2)
#' norm <- (norm-min(norm)) / max(norm-min(norm))
#' img <- outer(outer(norm, norm), norm)
#' img <- round(255 * img)
#' img[17:32,,] <- 255 - img[17:32,,]
#' img.nifti <- nifti(img) # create NIfTI object
#' 
#' writeNIfTI(img.nifti, "test-nifti-image-uint8", verbose=TRUE)
#' ## These files should be viewable in, for example, FSLview
#' ## Make sure you adjust the min/max values for proper visualization
#' data <- readNIfTI("test-nifti-image-uint8", verbose=TRUE)
#' image(img.nifti, oma=rep(2,4), bg="white")
#' image(data, oma=rep(2,4), bg="white")
#' abs.err <- abs(data - img.nifti)
#' image(as(abs.err, "nifti"), zlim=range(img.nifti), oma=rep(2,4),
#'       bg="white")
#' 
#' \dontrun{
#' ## Loop through all possible data types
#' datatypes <- list(code=c(2, 4, 8, 16, 64),
#'                   name=c("uint8", "int16", "int32", "float", "double"))
#' equal <- vector("list")
#' for (i in 1:length(datatypes$code)) {
#'   fname <- paste("test-nifti-image-", datatypes$name[i], sep="")
#'   rm(img.nifti)
#'   img.nifti <- nifti(img, datatype=datatypes$code[i])
#'   writeNIfTI(img.nifti, fname, verbose=TRUE)
#'   equal[[i]] <- all(readNIfTI(fname) == img)
#' }
#' names(equal) <- datatypes$name
#' unlist(equal)
#' }
#' @export
#' @rdname write_nifti
setMethod("writeNIfTI", signature(nim="nifti"), 
	  function(nim, filename, onefile=TRUE, gzipped=TRUE, verbose=FALSE,
                   warn=-1) {
            .writeNIfTI(nim, filename, onefile, gzipped, verbose, warn)
          })
#' @export
#' @rdname write_nifti
setMethod("writeNIfTI", signature(nim="anlz"), 
	  function(nim, filename, onefile=TRUE, gzipped=TRUE, verbose=FALSE,
                   warn=-1) {
            .writeNIfTI(as(nim, "nifti"), filename, onefile, gzipped,
                        verbose, warn)
          })
#' @export
#' @rdname write_nifti
setMethod("writeNIfTI", signature(nim="array"), 
	  function(nim, filename, onefile=TRUE, gzipped=TRUE, verbose=FALSE,
                   warn=-1) {
            .writeNIfTI(as(nim, "nifti"), filename, onefile, gzipped,
                        verbose, warn)
          })

.writeNIfTI <- function(nim, filename, onefile=TRUE, gzipped=TRUE,
                        verbose=FALSE, warn=-1) {
  ## Warnings?
  oldwarn <- getOption("warn")
  options(warn=warn)
  #### added so that range of the data will equal cal.min/cal.max
  nim <- calibrateImage(nim)
  ##### Added so that bad dimensions are dropped
  #   nim = drop_img_dim(nim)
  ## Basic error checking
  validNIfTI <- getValidity(getClassDef("nifti"))
  if (is.character(vnim <- validNIfTI(nim))) {
    stop(vnim)
  }
  ## Write header file...
  if (gzipped) {
    fid <- gzfile(paste(filename, "nii.gz", sep="."), "wb")
  } else {
    fid <- file(paste(filename, "nii", sep="."), "wb")
  }
  ## Extensions...
  extensions <- NULL
  if (is(nim, "niftiExtension")) {
    if (verbose) {
      cat("  niftiExtension detected!", fill=TRUE)
    }
    extensions <- nim@extensions
  }
  if (getOption("niftiAuditTrail") && is(nim, "niftiAuditTrail")) {
    if (verbose) {
      cat("  niftiAuditTrail detected!", fill=TRUE)
    }
    sec <- niftiAuditTrailToExtension(nim, getwd(), filename, match.call())
    extensions <- append(extensions, sec)
  }
  if (! is.null(extensions)) {
    ## update the vox_offset  FIXME twofile!
    totalesizes <- sum(unlist(lapply(extensions, function(x) x@"esize")))
    nim@"extender"[1] <- 1
    nim@"vox_offset" <- 352 + totalesizes
    if (verbose) {
      cat("  vox_offset =", nim@"vox_offset", fill=TRUE)
    }
  }
  ##
  writeBin(as.integer(nim@"sizeof_hdr"), fid, size=4)
  writeChar(nim@"data_type", fid, nchars=10, eos=NULL)
  writeChar(nim@"db_name", fid, nchars=18, eos=NULL)
  writeBin(as.integer(nim@"extents"), fid, size=4)
  writeBin(as.integer(nim@"session_error"), fid, size=2)
  writeChar(nim@"regular", fid, nchars=1, eos=NULL)
  writeBin(as.integer(nim@"dim_info"), fid, size=1)
  writeBin(as.integer(nim@"dim_"), fid, size=2)
  writeBin(as.double(nim@"intent_p1"), fid, size=4)
  writeBin(as.double(nim@"intent_p2"), fid, size=4)
  writeBin(as.double(nim@"intent_p3"), fid, size=4)
  writeBin(as.integer(nim@"intent_code"), fid, size=2)
  writeBin(as.integer(nim@"datatype"), fid, size=2)
  writeBin(as.integer(nim@"bitpix"), fid, size=2)
  writeBin(as.integer(nim@"slice_start"), fid, size=2)
  writeBin(as.double(nim@"pixdim"), fid, size=4)
  writeBin(as.double(nim@"vox_offset"), fid, size=4) # default offset = 352
  writeBin(as.double(nim@"scl_slope"), fid, size=4)
  writeBin(as.double(nim@"scl_inter"), fid, size=4)
  writeBin(as.integer(nim@"slice_end"), fid, size=2)
  writeBin(as.integer(nim@"slice_code"), fid, size=1)
  writeBin(as.integer(nim@"xyzt_units"), fid, size=1)
  writeBin(as.double(nim@"cal_max"), fid, size=4)
  writeBin(as.double(nim@"cal_min"), fid, size=4)
  writeBin(as.double(nim@"slice_duration"), fid, size=4)
  writeBin(as.double(nim@"toffset"), fid, size=4)
  writeBin(as.integer(nim@"glmax"), fid, size=4)
  writeBin(as.integer(nim@"glmin"), fid, size=4)
  writeChar(nim@"descrip", fid, nchars=80, eos=NULL)
  writeChar(nim@"aux_file", fid, nchars=24, eos=NULL)
  writeBin(as.integer(nim@"qform_code"), fid, size=2)
  writeBin(as.integer(nim@"sform_code"), fid, size=2)
  writeBin(as.double(nim@"quatern_b"), fid, size=4)
  writeBin(as.double(nim@"quatern_c"), fid, size=4)
  writeBin(as.double(nim@"quatern_d"), fid, size=4)
  writeBin(as.double(nim@"qoffset_x"), fid, size=4)
  writeBin(as.double(nim@"qoffset_y"), fid, size=4)
  writeBin(as.double(nim@"qoffset_z"), fid, size=4)
  writeBin(as.double(nim@"srow_x"), fid, size=4)
  writeBin(as.double(nim@"srow_y"), fid, size=4)
  writeBin(as.double(nim@"srow_z"), fid, size=4)
  writeChar(nim@"intent_name", fid, nchars=16, eos=NULL)
  writeChar(nim@"magic", fid, nchars=4, eos=NULL)
  writeBin(as.integer(nim@"extender"), fid, size=1)
  ## writeChar(as.character(nim@"extender"), fid, nchars=4, eos=NULL)
  ## Extensions?
  if (nim@"extender"[1] > 0 || nim@"vox_offset" > 352) {
    if (! is.null(extensions)) {
      if (verbose) {
        cat("  writing niftiExtension(s) at byte =", seek(fid), fill=TRUE)
      }
      lapply(extensions,
             function(x) {
               writeBin(as.integer(x@"esize"), fid, size=4)
               writeBin(as.integer(x@"ecode"), fid, size=4)
               ## Write out all the characters in the data section
	       writeChar(x@"edata", fid, nchars=nchar(x@"edata"), eos=NULL)
	       ## add margin to write \0 till 0 mod 16
	       margin <- (-(nchar(x@"edata", type="bytes") + 8) %% 16) 
	       if (margin > 0) {
		 writeBin(rep("", margin), fid, size=margin)
	       }
               invisible()
             })
    } else {
      stop("@extender set but", nim, "has no extensions.")
    }
  }
  ## reorient?
  if (nim@"reoriented") {
    data <- as.vector(inverseReorient(nim))
  } else {
    data <- as.vector(nim@.Data)
  }
  ## Write image file...
  if (verbose) {
    cat("  writing data at byte =", seek(fid), fill=TRUE)
  }
  switch(as.character(nim@"datatype"),
         "2" = writeBin(as.integer(data), fid, size=nim@"bitpix"/8),
         "4" = writeBin(as.integer(data), fid, size=nim@"bitpix"/8),
         "8" = writeBin(as.integer(data), fid, size=nim@"bitpix"/8),
         "16" = writeBin(as.double(data), fid, size=nim@"bitpix"/8),
         "64" = writeBin(as.double(data), fid, size=nim@"bitpix"/8),
         "512" = writeBin(as.integer(data), fid, size=nim@"bitpix"/8)
         )
  close(fid)
  ## Warnings?
  options(warn=oldwarn)
  invisible()
}

############################################################################
############################################################################
############################################################################
setGeneric("writeANALYZE", function(aim,  ...) standardGeneric("writeANALYZE"))
#' @title writeANALYZE
#' 
#' @description This function saves an Analyze-class object to a single binary file in
#' Analyze format.
#' 
#' @details The \code{writeANALYZE} function utilizes the internal \code{writeBin} and
#' \code{writeChar} command to write information to a binary file.
#' 
#' @name writeANALYZE-methods
#' @aliases writeANALYZE writeANALYZE-methods writeANALYZE,anlz-method
#' @docType methods
#' @param aim is an object of class \code{anlz}.
#' @param filename is the path and file name to save the Analyze file pair
#' (.hdr,img) \bold{without} the suffixes.
#' @param gzipped is a character string that enables exportation of compressed
#' (.gz) files (default = \code{TRUE}).
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @param warn is a number to regulate the display of warnings (default = -1).
#' See \code{\link{options}} for more details.
#' @return Nothing.
#' @section Methods: \describe{ \item{object = "anlz"}{Write ANALYZE volume to
#' disk.} }
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{writeAFNI}}, \code{\link{writeNIfTI}}
#' @references Analyze 7.5\cr \url{https://rportal.mayo.edu/bir/ANALYZE75.pdf}
#' @keywords file methods
#' @examples
#' 
#' norm <- dnorm(seq(-5, 5, length=32), sd=2)
#' norm <- (norm-min(norm)) / max(norm-min(norm))
#' img <- outer(outer(norm, norm), norm)
#' img <- round(255*img)
#' img[17:32,,] <- 255 - img[17:32,,]
#' img.anlz <- anlz(img) # create Analyze object
#' 
#' writeANALYZE(img.anlz, "test-anlz-image-uint8", verbose=TRUE)
#' ## These files should be viewable in, for example, FSLview
#' ## Make sure you adjust the min/max values for proper visualization
#' data <- readANALYZE("test-anlz-image-uint8", verbose=TRUE)
#' image(img.anlz, oma=rep(2,4), bg="white")
#' image(data, oma=rep(2,4), bg="white")
#' abs.err <- abs(data - img.anlz)
#' image(as(abs.err, "anlz"), zlim=range(img.anlz), oma=rep(2,4), bg="white")
#' 
#' \dontrun{
#' ## Loop through all possible data types
#' datatypes <- list(code=c(2, 4, 8, 16, 64),
#'                   name=c("uint8", "int16", "int32", "float", "double"))
#' equal <- vector("list")
#' for (i in 1:length(datatypes$code)) {
#'   fname <- paste("test-anlz-image-", datatypes$name[i], sep="")
#'   rm(img.anlz)
#'   img.anlz <- anlz(img, datatype=datatypes$code[i])
#'   writeANALYZE(img.anlz, fname)
#'   equal[[i]] <- all(readANALYZE(fname) == img)
#' }
#' names(equal) <- datatypes$name
#' unlist(equal)
#' }
#' @export
#' @rdname write_anlz 
setMethod("writeANALYZE", signature(aim="anlz"), 
	  function(aim, filename, gzipped=TRUE, verbose=FALSE, warn=-1) {
            .writeANALYZE(aim, filename, gzipped, verbose, warn)
          })
.writeANALYZE <- function(aim, filename, gzipped=TRUE, verbose=FALSE,
                          warn=-1) {
  ## Warnings?
  oldwarn <- getOption("warn")
  options(warn=warn)
  ## Basic error checking
  validANALYZE <- getValidity(getClassDef("anlz"))
  if (is.character(error <- validANALYZE(aim))) {
    cat(error, fill=TRUE)
    stop("-- aborting writeANALYZE --")
  }
  ## Write header file...
  if (gzipped) {
    fid <- gzfile(paste(filename, ".hdr.gz", sep=""), "wb")
  } else {
    fid <- file(paste(filename, ".hdr", sep=""), "wb")
  }
  ## header_key
  writeBin(as.integer(aim@"sizeof_hdr"), fid, size=4)  #  0 + 4
  writeChar(aim@"data_type", fid, nchars=10, eos=NULL) #  4 + 10
  writeChar(aim@"db_name", fid, nchars=18, eos=NULL)   #  14 + 18
  writeBin(as.integer(aim@"extents"), fid, size=4)     #  32 + 4
  writeBin(as.integer(aim@"session_error"), fid, size=2) # 36 + 2
  writeChar(aim@"regular", fid, nchars=1, eos=NULL)    # 38 + 1
  writeChar(aim@"hkey_un0", fid, nchars=1, eos=NULL)   # 39 + 1
  ##                                                   40 bytes
  ## image_dimension
  writeBin(as.integer(aim@"dim_"), fid, size=2)        # 0 + (2 x 8)
  writeChar(aim@"vox_units", fid, nchars=4, eos=NULL)  # 16 + 4
  writeChar(aim@"cal_units", fid, nchars=8, eos=NULL)  # 20 + 8
  writeBin(as.integer(aim@"unused1"), fid, size=2)     # 28 + 2
  writeBin(as.integer(aim@"datatype"), fid, size=2)    # 30 + 2
  writeBin(as.integer(aim@"bitpix"), fid, size=2)      # 32 + 2
  writeBin(as.integer(aim@"dim_un0"), fid, size=2)     # 34 + 2
  writeBin(aim@"pixdim", fid, size=4)                  # 36 + (4 x 8)
  writeBin(aim@"vox_offset", fid, size=4)              # 68 + 4
  writeBin(aim@"funused1", fid, size=4)                # 72 + 4
  writeBin(aim@"funused2", fid, size=4)                # 76 + 4
  writeBin(aim@"funused3", fid, size=4)                # 80 + 4
  writeBin(as.double(aim@"cal_max"), fid, size=4)     # 84 + 4
  writeBin(as.double(aim@"cal_min"), fid, size=4)     # 88 + 4
  writeBin(as.integer(aim@"compressed"), fid, size=4)  # 92 + 4
  writeBin(as.integer(aim@"verified"), fid, size=4)    # 96 + 4
  writeBin(as.integer(aim@"glmax"), fid, size=4)       # 100 + 4
  writeBin(as.integer(aim@"glmin"), fid, size=4)       # 104 + 4
  ##                                                   108 bytes
  ## data_history
  writeChar(aim@"descrip", fid, nchars=80, eos=NULL)   # 0 + 80
  writeChar(aim@"aux_file", fid, nchars=24, eos=NULL)  # 80 + 24
  writeChar(aim@"orient", fid, nchars=1, eos=NULL)     # 104 + 1
  writeBin(as.integer(aim@"origin"), fid, size=2)      # 105 + (2 x 5)
  writeChar(aim@"generated", fid, nchars=10, eos=NULL) # 115 + 10
  writeChar(aim@"scannum", fid, nchars=10, eos=NULL)   # 125 + 10
  writeChar(aim@"patient_id", fid, nchars=10, eos=NULL) # 135 + 10
  writeChar(aim@"exp_date", fid, nchars=10, eos=NULL)  # 145 + 10
  writeChar(aim@"exp_time", fid, nchars=10, eos=NULL)  # 155 + 10
  writeChar(aim@"hist_un0", fid, nchars=3, eos=NULL)   # 165 + 3
  writeBin(as.integer(aim@"views"), fid, size=4)       # 168 + 4
  writeBin(as.integer(aim@"vols_added"), fid, size=4)  # 172 + 4
  writeBin(as.integer(aim@"start_field"), fid, size=4) # 176 + 4
  writeBin(as.integer(aim@"field_skip"), fid, size=4)  # 180 + 4
  writeBin(as.integer(aim@"omax"), fid, size=4)        # 184 + 4
  writeBin(as.integer(aim@"omin"), fid, size=4)        # 188 + 4
  writeBin(as.integer(aim@"smax"), fid, size=4)        # 192 + 4
  writeBin(as.integer(aim@"smin"), fid, size=4)        # 196 + 4
  ##                                                   200 bytes
  close(fid)
  ## Write image file...
  if (gzipped) {
    fid <- gzfile(paste(filename, "img.gz", sep="."), "wb")
  } else {
    fid <- file(paste(filename, "img", sep="."), "wb")
  }
  dims <- 2:(1+aim@"dim_"[1])
  if (verbose) {
    cat("  dims =", aim@"dim_"[dims], fill=TRUE)
  }
  data <- as.vector(aim@.Data)
  switch(as.character(aim@"datatype"),
         "1" = writeBin(as.integer(data), fid, size=aim@"bitpix"/8),
         "2" = writeBin(as.integer(data), fid, size=aim@"bitpix"/8),
         "4" = writeBin(as.integer(data), fid, size=aim@"bitpix"/8),
         "8" = writeBin(as.integer(data), fid, size=aim@"bitpix"/8),
         "16" = writeBin(as.numeric(data), fid, size=aim@"bitpix"/8),
         "64" = writeBin(as.double(data), fid, size=aim@"bitpix"/8))
  close(fid)
  ## Warnings?
  options(warn=oldwarn)
  invisible()
}
