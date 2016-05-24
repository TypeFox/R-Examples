##
## Copyright (c) 2010-2011 Brandon Whitcher
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

#' Construct Data Frame from DICOM Headers
#' 
#' A data frame is created given the valid DICOM fields provided by the user.
#' 
#' 
#' @param hdrs List object of DICOM headers.
#' @param stringsAsFactors Logical variable to be passed to \code{data.frame}.
#' @param collapse Character string used to \code{paste} DICOM group, element
#' and value fields.
#' @param colSort Logical variable (default = \code{TRUE}) to sort column names
#' in the table.
#' @param verbose Flag to provide text-based progress bar (default =
#' \code{FALSE}).
#' @param debug Logical variable (default = \code{FALSE}) that regulates to
#' display of intermediate processing steps.
#' @return Data frame where the rows correspond to images and the columns
#' correspond to the UNION of all DICOM fields across all files in the list.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @references Whitcher, B., V. J. Schmid and A. Thornton (2011).  Working with
#' the DICOM and NIfTI Data Standards in R, \emph{Journal of Statistical
#' Software}, \bold{44} (6), 1--28.  \url{http://www.jstatsoft.org/v44/i06}
#' 
#' Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}
#' @keywords misc
#' @export dicomTable
dicomTable <- function(hdrs, stringsAsFactors=FALSE, collapse="-",
                       colSort=TRUE, verbose=FALSE, debug=FALSE) {
  myMerge <- function(df1, df2) {
    if (anyDuplicated(names(df1)) != 0) {
      warning("Duplicated group-element tags have been removed!")
      df1 <- df1[, ! duplicated(names(df1))]
    }
    if (anyDuplicated(names(df2)) != 0) {
      warning("Duplicated group-element tags have been removed!")
      df2 <- df2[, ! duplicated(names(df2))]
    }
    if (! all(names(df2) %in% names(df1))) {
      newCols <- names(df2)[! names(df2) %in% names(df1)]
      ## newcols <- setdiff(names(df2), names(df1)) # removes duplicates!
      newDf <- as.data.frame(lapply(newCols, function(i, x) rep(NA, x),
                                    x = nrow(df1)))
      names(newDf) <- newCols
      df1 <- cbind(df1, newDf)
    }
    if (! all(names(df1) %in% names(df2))) {
      newCols <- names(df1)[! names(df1) %in% names(df2)]
      ## newCols <- setdiff(names(df1), names(df2)) # removes duplicates!
      newDf <- as.data.frame(lapply(newCols, function(i, x) rep(NA, x),
                                    x = nrow(df2)))
      names(newDf) <- newCols
      df2 <- cbind(df2, newDf)
    }
    rbind(df1, df2)
  }
  ## Use first record to establish data.frame
  csv <- data.frame(matrix(hdrs[[1]]$value, 1, nrow(hdrs[[1]])),
                    stringsAsFactors=stringsAsFactors)
  names(csv) <-
    paste(sub("^-", "", gsub("[^0-9]+", "-", hdrs[[1]]$sequence)),
          as.vector(apply(hdrs[[1]][,1:3], 1, paste, collapse=collapse)),
          sep="")
  ## Loop through all records and "merge" them
  if ((nhdrs <- length(hdrs)) > 1) {
    if (verbose) {
      cat(" ", nhdrs, "files to be processed by dicomTable()", fill=TRUE)
      tpb <- txtProgressBar(min=0, max=nhdrs, style=3)
    }
    for (l in 2:nhdrs) {
      if (debug) {
        cat("  l =", l, fill=TRUE)
      }
      if (verbose) {
        setTxtProgressBar(tpb, l)
      }
      temp <- data.frame(matrix(hdrs[[l]]$value, 1, nrow(hdrs[[l]])),
                         stringsAsFactors=stringsAsFactors)
      names(temp) <-
        paste(sub("^-", "", gsub("[^0-9]+", "-", hdrs[[l]]$sequence)),
              as.vector(apply(hdrs[[l]][,1:3], 1, paste, collapse=collapse)),
              sep="")
      old.nrow <- nrow(csv)
      csv <- myMerge(csv, temp)
      if (nrow(csv) == old.nrow) {
        warning("Duplicate row was _not_ inserted in data.frame (csv)")
        csv <- rbind(csv, NA)
      }
    }
    if (verbose) {
      close(tpb)
    }
    row.names(csv) <- names(hdrs)
  }
  if (colSort) {
    return(csv[, order(names(csv))])
  } else {
    return(csv)
  }
}

#' Extract Single Field from DICOM Headers
#' 
#' A particular DICOM field is extracted for a collection of DICOM headers.
#' 
#' The DICOM field is extracted from each DICOM header and placed into a
#' vector.
#' 
#' @param hdrs List object of DICOM headers.
#' @param string DICOM field name.
#' @param numeric Logical; values are converted to numbers when \code{TRUE}.
#' @param names Logical; file names are kept with elements of the vector.
#' @param inSequence Logical; whether or not to look into SequenceItem
#' elements.
#' @return Vector of values from the requested DICOM field.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{readDICOM}}
#' @references Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}
#' @keywords misc
#' @examples
#' 
#' x <- readDICOMFile(system.file("dcm/Abdo.dcm", package="oro.dicom"))
#' seriesDescription <- extractHeader(x$hdr, "SeriesDescription", numeric=FALSE)
#' IOP <- extractHeader(x$hdr, "ImageOrientationPatient", numeric=FALSE)
#' 
#' @export extractHeader
extractHeader <- function(hdrs, string, numeric=TRUE, names=FALSE,
                          inSequence=TRUE) {
  if (is.data.frame(hdrs)) {
    hdrs <- list(hdrs)
  }
  out.list <- lapply(hdrs,
                     function(hdr, string, inSequence) {
                       if (inSequence) {
                         sequence <- FALSE
                       } else {
                         sequence <- nchar(hdr$sequence) > 0
                       }
                       index <- which(hdr$name %in% string & !sequence)
                       if (sum(index) > 0) {
                         hdr$value[index]
                       } else {
                         NA
                       }
                     }, string=string, inSequence=inSequence)
  out.names <- names(out.list)
  out.vec <- unlist(out.list)
  if (numeric) {
    out.vec <- as.numeric(out.vec)
  }
  if (names) {
    names(out.vec) <- out.names
  } else {
    out.vec <- as.vector(out.vec)
  }
  return(out.vec)
}

#' Converts DICOM Header Field to a Matrix
#' 
#' Converts a vector of DICOM header information, assuming there are multiple
#' entries per element of the vector, into a matrix.
#' 
#' 
#' @param hdr is the result from extracting information from a DICOM header
#' field; e.g., using \code{extractHeader}.
#' @param ncol is the number of columns.
#' @param sep is the character string required to split entries in the header
#' field.
#' @param byrow is a logical variable (default = \code{TRUE}) telling the
#' routine to populate the matrix by rows then columns.
#' @return Matrix with \code{length(hdr)} rows and \code{ncol} columns.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{extractHeader}}, \code{\link{matrix}}
#' @references Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}
#' @keywords misc
#' @examples
#' 
#' x <- readDICOMFile(system.file("dcm/Abdo.dcm", package="oro.dicom"))
#' pixelSpacing <- extractHeader(x$hdr, "PixelSpacing", numeric=FALSE)
#' pSmat <- header2matrix(pixelSpacing, ncol=2)
#' IOP <- extractHeader(x$hdr, "ImageOrientationPatient", numeric=FALSE)
#' IOPmat <- header2matrix(IOP, ncol=6)
#' 
#' @export header2matrix
header2matrix <- function(hdr, ncol, sep=" ", byrow=TRUE) {
  matrix(as.numeric(unlist(strsplit(hdr, sep))), ncol=ncol, byrow=byrow)
}



#' Match String to DICOM Header Field
#' 
#' A convenient wrapper function that utilizes internal functions to match
#' character strings with the DICOM header information.
#' 
#' 
#' @param hdr is the result from extracting information from a DICOM header
#' field; e.g., using \code{extractHeader}.
#' @param string is a character string to be matched with the DICOM header.
#' @return A logical vector of length \code{length(hdr)}.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{extractHeader}}
#' @references Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}
#' @examples
#' 
#' x <- readDICOMFile(system.file("dcm/Abdo.dcm", package="oro.dicom"))
#' modality <- extractHeader(x$hdr, "Modality", numeric=FALSE)
#' matchHeader(modality, "mr") # case insensitive by default
#' 
#' @export matchHeader
matchHeader <- function(hdr, string) {
  ifelse(is.na(hdr), FALSE, regexpr(string, hdr, ignore.case=TRUE) > -1)
}

#' Write DICOM Table to ASCII File
#' 
#' A wrapper to \code{write.table} specifically for DICOM tables.
#' 
#' This function is a straightforward wrapper to \code{write.table}.
#' 
#' @param dtable The DICOM table.
#' @param filename Name of the file to be created.
#' @param ... Additional parameters to be passed to \code{write.table}.
#' @return None.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{write.table}}
#' @references Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}
#' @keywords file
#' @export writeHeader
writeHeader <- function(dtable, filename, ...) {
  write.table(dtable, filename, quote=FALSE, sep="\t", ...)
}

#' Check String Against DICOM Header Field to Produce Error Message or NEXT
#' 
#' A function designed to \code{break} out of loops given information (or the
#' lackthereof) contained in the DICOM header.
#' 
#' 
#' @param dcm is the DICOM list structure.
#' @param string is a character string to be matched with the DICOM header.
#' @param reference is the scalar/vector of character strings to check against
#' the DICOM header output.
#' @param str.warning is a text string for the warning.
#' @param htmlfile is the \pkg{hwriter} object for the HTML file (default =
#' \code{NULL}.
#' @param heading is the HTML tag <H?> (default = \code{3}).
#' @param numeric is the argument to be passed to \code{matchHeader}.
#' @return An expression to be evaluated and HTML content.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{extractHeader}}, \code{\link{matchHeader}}
#' @references Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}
#' @export nextHeader
nextHeader <- function(dcm, string, reference, str.warning,
                       htmlfile=NULL, heading=3, numeric=FALSE) {
  header <- extractHeader(dcm$hdr, string=string, numeric=numeric)
  for (i in 1:length(reference)) {
    if (any(matchHeader(header, string=reference[i]))) {
      if (! is.null(htmlfile)) {
        requireNamespace("hwriter", quietly=TRUE)
        hwriter::hwrite(str.warning, htmlfile, heading=3)
      } else {
        warning(str.warning)
      }
      return(expression(next))
    }
  }
  invisible()
}
