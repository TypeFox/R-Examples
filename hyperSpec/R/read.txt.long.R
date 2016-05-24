###-----------------------------------------------------------------------------
###
###  read.txt.long: import measurements from .txt file
###
###  Format:
###  (y x) wl int
###

##' Import and Export of hyperSpec objects
##' Besides \code{\link[base]{save}} and \code{\link[base]{load}}, two general
##' ways to import and export data into \code{hyperSpec} objects exist.
##' 
##' Firstly, hyperSpec objects can be imported and exported as ASCII files.
##' 
##' A second option is using the package \code{\link[R.matlab]{R.matlab}} which
##' provides the functions \code{\link[R.matlab]{readMat}} and
##' \code{\link[R.matlab]{writeMat}}.
##' 
##' hyperSpec comes with a number of pre-defined functions to import
##' manufacturer specific file formats. For details, see \code{vignette
##' ("file-io")}.
##' 
##' \code{\link[hyperSpec]{read.spc}} imports Thermo Galactic's .spc file
##' format, and ENVI files may be read using
##' \code{\link[hyperSpec]{read.ENVI}}.
##' 
##' These functions are very flexible and provide lots of arguments.
##' 
##' If you use them to read or write manufacturer specific ASCII formats,
##' please consider writing a wrapper function and contributing this function
##' to \pkg{hyperSpec}.  An example is in the \dQuote{flu} vignette (see
##' \code{vignette ("flu", package = "hyperSpec"}).
##' 
##' Note that R accepts many packed formats for ASCII files, see
##' \code{\link[base]{connections}}. For .zip files, see
##' \code{\link[utils]{unzip}}.
##' 
##' For further information, see the examples below and the documentation of
##' \code{\link[R.matlab]{R.matlab}}.
##' 
##' @aliases read.txt.long import export
##' @param file filename or connection
##' @param cols the column names specifying the column order.
##' 
##' For data import, a list with elements \code{colname = label}; for export a
##'   character vector with the colnames.  Use \code{wavelength} to specify the
##'   wavelengths.
##' @param header the file has (shall have) a header line
##' @param ... arguments handed to \code{\link[utils]{read.table}} and
##'   \code{\link[utils]{write.table}}, respectively.
##' @param decreasing logical vector giving the sort order
##' @author C. Beleites
##' @seealso \code{\link[utils]{read.table}} and
##'   \code{\link[utils]{write.table}}
##' 
##' \code{\link[R.matlab]{R.matlab}} for .mat files
##' 
##' \code{\link[hyperSpec]{read.ENVI}} for ENVI data
##' 
##' \code{\link[hyperSpec]{read.spc}} for .spc files
##' 
##' Manufacturer specific file formats: \code{\link{scan.txt.Renishaw}}
##' @rdname textio
##' @keywords IO file
##' @export
##' @examples
##' 
##' 
##' \dontrun{vignette  ("file-io")}
##' 
##' ## export & import matlab files
##' if (require (R.matlab)){
##'    # export to matlab file
##'    writeMat ("test.mat", x = flu[[]], wavelength = flu@@wavelength,
##'              label = lapply (flu@@label, as.character))
##' 
##'    # reading a matlab file
##'    data <- readMat ("test.mat")
##'    print (data)
##'    mat <- new ("hyperSpec", spc = data$x,
##'                wavelength = as.numeric(data$wavelength),
##'                label = data$label[,,1])
##' }
##' 
##' ## ascii export & import
##' 
##' 
##' write.txt.long (flu,  file = "flu.txt", cols = c(".wavelength", "spc", "c"), 
##' 		order = c("c", ".wavelength"),
##' 		decreasing = c(FALSE, TRUE))
##' 
##' read.txt.long (file = "flu.txt", cols = list (.wavelength = expression (lambda / nm), 
##'       spc= "I / a.u", c = expression ("/" (c, (mg/l)))))
##' 
##' write.txt.wide (flu,  file = "flu.txt", cols = c("c", "spc"), 
##' 		col.labels = TRUE, header.lines = 2, row.names = TRUE)
##' 
##' write.txt.wide (flu,  file = "flu.txt", col.labels = FALSE, row.names = FALSE)
##' 
##' read.txt.wide (file = "flu.txt", 
##'       cols = list (c=expression ("/"("c", "mg/l")), spc="I / a.u", .wavelength = "lambda / nm"),
##' 		header = TRUE)
##' 
##' 
read.txt.long <- function (file = stop ("file is required"),
                           cols = list (
                             .wavelength = expression (lambda / nm),
                             spc = "I / a.u."),
                           header = TRUE,
                           ...){
  txtfile <- read.table (file = file, header = header, ...)

  if (header){
    cln <- match (colnames (txtfile), names (cols))
    cln <- cols[cln]
    names (cln) <- colnames (txtfile)
    cols <- cln
    rm (cln)
  } else {
    if (ncol (txtfile) != length (cols)){
      warning (paste ("cols does not correspond to the columns in", file,
                      ". Guessing remaining columns."))
      cols <- c (character (ncol (txtfile) - 2), cols)
    }
  }


  if (is.na (match ("spc", names (cols))))
    stop ("cols$spc must exist.")

  wavelength <- match (".wavelength", names (cols))
  if (is.na (wavelength))
    stop ("cols$.wavelength must exist.")

  colnames (txtfile) <- names (cols)

  ## wavelength axis
  wavelength <- as.numeric (levels (as.factor (txtfile$.wavelength)))
  
  spc <- as.matrix (unstack (txtfile, form = spc ~ .wavelength))
  if ((nrow (spc)  == length (wavelength)) & (ncol (spc) != length (wavelength)))
    spc <- t (spc)

  colnames (spc) <- levels (txtfile$.wavelength)

  txtfile <- txtfile [txtfile$.wavelength == txtfile$.wavelength[1], ]
  txtfile$.wavelength <- NULL
  txtfile$spc <- I (spc)

  new ("hyperSpec", wavelength = wavelength, data = txtfile, labels = cols)
}
