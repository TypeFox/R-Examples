##' Import Raman Spectra/Maps from Witec Instrument via ASCII files
##'
##' \code{scan.dat.Witec} reads Witec's ASCII exported data which comes in separate files with x and
##' y data. \code{scan.txt.Witec} reads Witec ASCII files where the first column gives the wavelength
##' axes and the other columns the spectra.
##' @title File Import Witec Raman
##' @param filex filename wavelength axis file
##' @param filey filename intensity file
##' @param points.per.line number of spectra in x direction of the map
##' @param lines.per.image number of spectra in y direction
##' @param ...,quiet handed to \code{\link[base]{scan}}
##' @param nwl number of wavelengths, if \code{NULL}, \code{readLines} is used to determine
##' \code{nwl} automatically.
##' @param remove.zerospc WiTEC Control saves spectra consisting of zeros only if e.g. a map was
##' aborted. The default is to remove these spectra.
##' @return a hyperSpec object
##' @author Claudia Beleites
##' @seealso \code{vignette ("fileio")} for more information on file import and
##' 
##' \code{\link{options}} for details on options. 
##' @export 
scan.dat.Witec <- function (filex = stop ("filename or connection needed"),
                            filey = sub ("-x", "-y", filex),
                            points.per.line = NULL,
                            lines.per.image = NULL,
                            ..., 
                            quiet = hy.getOption ("debuglevel") < 1L){
  wl <- scan (file = filex, ..., quiet = quiet)
  spc <- scan (file = filey, ..., quiet = quiet)

  dim (spc) <- c (length (wl), length (spc) / length (wl))

  spc <- new ("hyperSpec", wavelength = wl, spc = t (spc))

  if (!is.null (points.per.line))
    spc@data$x <- rep (seq_len (points.per.line), lines.per.image)

  if (!is.null (points.per.line))
    spc@data$y <- rep (- seq_len (lines.per.image), each = points.per.line)
    
  ## consistent file import behaviour across import functions
  .fileio.optional (spc, filey)
}

##' @rdname scan.dat.Witec
##' @param file filename or connection to ASCII file
##' @export
scan.txt.Witec <- function (file = stop ("filename or connection needed"),
                            points.per.line = NULL,
                            lines.per.image = NULL,
                            nwl = 1024,
                            remove.zerospc = TRUE,
                            ...){

  if (is.null (nwl)){
    txt <- readLines (file)
    nwl <- length (txt)
    txt <- scan (text = txt, ...)
  } else {
    txt <- scan (file, ...)
  }

  dim (txt) <- c (length (txt) / nwl, nwl)

  spc <- new ("hyperSpec", wavelength = txt [1, ], spc = txt [-1, ])

  if (!is.null (points.per.line))
    spc@data$x <- rep (seq_len (points.per.line), lines.per.image)

  if (!is.null (points.per.line))
    spc@data$y <- rep (- seq_len (lines.per.image), each = points.per.line)

  ## consistent file import behaviour across import functions
  .fileio.optional (spc, file)
}
