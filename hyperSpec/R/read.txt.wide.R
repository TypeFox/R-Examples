###-----------------------------------------------------------------------------
###
###  read.txt.wide
###
###  Format:
###  x y ... int (wl1)  int (wl2) ... int (wl p) z ...
###
##' Import/export of hyperSpec objects to/from ASCII files
##' A detailed discussion of hyperSpec's file import and export capabilities is given in vignette \quote{fileio}.
##'
##' Besides \code{\link[base]{save}} and \code{\link[base]{load}}, two general ways to import and
##' export data into \code{hyperSpec} objects exist.
##' 
##' Firstly, hyperSpec objects can be imported  and exported as ASCII files.
##' 
##'   A second option is using the package \code{\link[R.matlab]{R.matlab}}
##'   which provides the functions \code{\link[R.matlab]{readMat}} and
##'   \code{\link[R.matlab]{writeMat}}.
##' 
##'   hyperSpec comes with a number of pre-defined functions to import
##'   manufacturer specific file formats. For details, see \code{vignette
##'   ("fileio")}.
##'   
##'   \code{\link[hyperSpec]{read.spc}} imports Thermo Galactic's .spc file
##'   format, and ENVI files may be read using
##'   \code{\link[hyperSpec]{read.ENVI}}.  
##' 
##' These functions are very flexible and provide lots of arguments. 
##' 
##' If you use them to read or write manufacturer specific ASCII formats,
##' please consider writing a wrapper function and contributing this
##' function to \pkg{hyperSpec}.  An example is in the \dQuote{flu} vignette
##' (see \code{vignette ("flu", package = "hyperSpec"}).
##' 
##' Note that R accepts many packed formats for ASCII files, see
##' \code{\link[base]{connections}}. For .zip files, see \code{\link[utils]{unzip}}.
##' 
##' For further information, see the examples below, \code{vignette ("fileio")} and the documentation
##' of \code{\link[R.matlab]{R.matlab}}.
##' @seealso \code{vignette ("fileio")} and \url{http://hyperspec.r-forge.r-project.org/blob/fileio.pdf},
##' respectively
##' @aliases read.txt.wide 
##' @rdname textio
##' @param check.names handed to \code{\link[utils]{read.table}}. Make sure this is \code{FALSE}, if
##' the column names of the spectra are the wavelength values.
##' @export
read.txt.wide <- function (file = stop ("file is required"),
                           cols = list (
                             spc = "I / a.u.",
                             .wavelength = expression (lambda / nm)),
                           sep = '\t',
                           row.names = NULL,
                           check.names = FALSE,
                           ...){
  
  .wavelength <- match (".wavelength", names (cols))
  if (is.na (.wavelength))
    cols <- as.list (c (cols, .wavelength = expression (lambda / nm)))
  else
    if (.wavelength != length (cols))   # .wavelength should be at the end of cols
      cols <- cols [c (seq_along (cols)[-.wavelength], .wavelength)]

  ## columns containing the spectra
  spc <- match ("spc", names (cols))
  if (is.na (spc))
    stop ("cols$spc must exist.")

  txtfile <- read.table (file = file, check.names = check.names, row.names = row.names,
                         sep = sep, ...)

  spc <- 0 : (ncol (txtfile) - length (cols) + 1) + spc

  spc.data <- as.matrix (txtfile[, spc])
  txtfile <- txtfile [, -spc, drop = FALSE]
  txtfile$spc <- I (spc.data)

  ## enforce colnames given by cols
  colnames (txtfile) <- head (names (cols), -1)

  new ("hyperSpec", data = txtfile, labels = cols)
}
