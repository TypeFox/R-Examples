##' import Raman measurements from Renishaw ASCII-files
##' import Raman measurements from Renishaw (possibly compressed) .txt file.
##' 
##' The file may be of any file type that can be read by
##' \code{\link[base]{gzfile}} (i.e. text, or zipped by gzip, bzip2, xz or
##' lzma). .zip zipped files need to be read using \code{scan.zip.Renishaw}.
##' 
##' Renishaw .wxd files are converted to .txt ASCII files by their batch
##' converter. They come in a "long" format with columns (y x | time | z)?
##' wavelength intensity.  The first columns depend on the data type.
##' 
##' The corresponding possibilities for the \code{data} argument are:
##' \tabular{lll}{ \code{data} \tab columns \tab \cr \code{"spc"} \tab wl int
##' \tab single spectrum \cr \code{"zspc"}, \code{"depth"} \tab z wl int \tab
##' depth profile\cr \code{"ts"} \tab t wl int \tab time series\cr
##' \code{"xyspc"} \tab y x wl int \tab 2d map\cr }
##' 
##' This function allows reading very large ASCII files, but it does not work
##' on files with missing values (\code{NA}s are allowed).
##' 
##' If the file is so large that it sould be read in chunks and \code{nspc} is
##' not given, \code{scan.txt.Renishaw} tries to guess it by using \code{wc}
##' (if installed).
##' 
##' @aliases scan.txt.Renishaw scan.zip.Renishaw
##' @param file file name or connection
##' @param data type of file, one of "spc", "xyspc", "zspc", "depth", "ts", see
##'   details.
##' @param nlines number of lines to read in each chunk, if 0 or less read
##'   whole file at once.
##' 
##' \code{nlines} must cover at least one complete spectrum,i.e. \code{nlines}
##'   must be at least the number of data points per spectrum. Reasonable
##'   values start at \code{1e6}.
##' @param nspc number of spectra in the file
##' @param ... Arguments for \code{scan.txt.Renishaw}
##' @return the \code{hyperSpec} object
##' @rdname scan-txt-Renishaw
##' @export
##' @author C. Beleites
##' @seealso \code{\link{read.txt.long}}, \code{\link{read.txt.wide}},
##'   \code{\link[base]{scan}}
##' @keywords IO file
scan.txt.Renishaw <- function (file = stop ("file is required"),
                               data = "xyspc", nlines = 0, nspc = NULL){
  cols <- switch (data,
                  spc = NULL,
                  xyspc = list (y = expression ("/" (y, mu * m)), 
                    x = expression ("/" (x, mu * m))), 
                  zspc = ,
                  depth = list (z = expression ("/" (z, mu * m))),
                  ts = 	list (t = "t / s"),
                  stop ("unknown format for Renishaw .txt files.")
                  )
  cols <- c  (cols, list (.wavelength = expression (Delta * tilde(nu) / cm^-1) ,
                          spc = "I / a.u."))

  if (!is (file, "connection"))
    file <- gzfile (file, "r")

  on.exit(close(file))

  first <- scan(file, nlines = 1, quiet = TRUE)
  
  ncol <- length (first)

  if (ncol == 0)
    return (new ("hyperSpec"))

  if (ncol != length (cols))
    stop (paste ("File has", ncol, "columns, while 'cols' gives", length (cols)))

  fbuf <- matrix (c (first, scan (file, quiet = TRUE, nlines = nlines)),
                     ncol = ncol,  byrow = TRUE)

  ## wavelength axis
  wl <- rep (TRUE,  nrow (fbuf))
  for (i in seq_len (ncol (fbuf) - 2))
    wl [wl] <- fbuf [wl, i] == fbuf [1, i]

  wl <- fbuf[wl, ncol - 1]

  ## if the file is to be read in chunks
  ## try to find out how many lines it has
  if (is.null (nspc))
    if (nlines > 0){
      nspc <- wc (summary(file)$description, "lines")
      if (is.null (nspc))
        stop ("failed guessing nspc.")
      else {
        cat ("Counted", nspc[1,1], "lines or ")
        nspc <- nspc[1,1] / length (wl)
        cat (nspc, "spectra.\n")
      }
    } else {
      nspc <- nrow (fbuf) / length (wl)
    }

  data <- matrix (NA, ncol = ncol - 2, nrow = nspc)
  colnames (data) <- head (names (cols), -2)
  pos.data <- 0

  spc <- numeric (nspc * length (wl))
  pos.spc <- 0

  while (length (fbuf > 0)){
    if (nlines > 0) cat (".")
    spc [pos.spc + seq_len (nrow (fbuf))] <- fbuf [, ncol]
    pos.spc <- pos.spc + nrow (fbuf)

    tmp <- fbuf [fbuf[, ncol - 1] == wl [1], seq_len (ncol - 2), drop = FALSE]

    data [pos.data + seq_len (nrow (tmp)), ] <- tmp
    pos.data <- pos.data + nrow (tmp)

    fbuf <- matrix (scan (file, quiet = TRUE, nlines = nlines), ncol = ncol,
                    byrow = TRUE)

    if (length (fbuf > 0) & ! all(unique (fbuf[, ncol - 1]) %in% wl))
      stop ("Wavelengths do not correspond to that of the other chunks.",
            "Is the size of the first chunk large enough to cover a complete",
            "spectrum?")
  }
  if (nlines > 0) cat ("\n")

  spc <- matrix (spc, ncol = length (wl), nrow = nspc, byrow = TRUE)

  orderwl (new ("hyperSpec", spc = spc, data = as.data.frame (data),
                wavelength = wl, label = cols))
}

##' @export
##' @param txt.file name of the .txt file in the .zip archive. Defaults to zip
##'   file's name with suffix .txt instead of .zip
##' @rdname scan-txt-Renishaw
scan.zip.Renishaw <- function (file = stop ("filename is required"),
                               txt.file = sub ("[.]zip", ".txt", basename (file)), ...)
  scan.txt.Renishaw (file = unz (file, filename = txt.file, "r"), ...)
