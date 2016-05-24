##' JCAMP-DX Import for Shimadzu Library Spectra
##'
##' this is a first rough import function for JCAMP-DX spectra.
##'
##' So far, AFFN and PAC formats are supported for simple XYDATA, DATA TABLEs and PEAK TABLEs.
##'
##' NTUPLES / PAGES are not (yet) supported.
##'
##' DIF, DUF, DIFDUP and SQZ data formats are not (yet) supported. 
##' 
##' @note JCAMP-DX support is incomplete and the functions may change without notice. See
##' \code{vignette ("fileio")}  and the details section.
##' @param filename file name and path of the .jdx file
##' @param encoding encoding of the JCAMP-DX file (used by \code{\link[base]{readLines}})
##' @param header list with manually set header values
##' @param keys.hdr2data index vector indicating which header entries should be tranfered into the
##' extra data. Usually a character vector of labels (lowercase, without and dashes, blanks,
##' underscores). If \code{TRUE}, all header entries are read.
##' @param ... further parameters handed to the data import function, e.g.
##' \tabular{ll}{
##' \code{xtol} \tab tolerance for checking calculated x values against checkpoints at beginning
##'                  of line, defaults to XFACTOR\cr
##' \code{ytol} \tab tolerance for checking Y values against MINY and MAXY, defaults to YFACTOR\cr
##' }
##' @param NA.symbols character vector of text values that should be converted to \code{NA}
##' @param collapse.multi should hyperSpec objects from multispectra files be collapsed into one
##' hyperSpec object (if \code{FALSE}, a list of hyperSpec objects is returned).
##' @return hyperSpec object
##' @author C. Beleites with contributions by Bryan Hanson
##' @export
read.jdx <- function(filename = stop ("filename is needed"), encoding = "",
                     header = list (), keys.hdr2data = FALSE, ...,
                     NA.symbols = c ("NA", "N/A", "N.A."),
                     collapse.multi = TRUE){
  
  ## see readLines help: this way, encoding is translated to standard encoding on current system.
  file <- file (filename, "r", encoding = encoding, blocking = FALSE)
  jdx <- readLines (file)
  close (file)
  
  ## start & end of spectra header and data
  hdrstart <- grep ("^[[:blank:]]*##TITLE=", jdx)
  if (length (hdrstart) == 0L) stop ("No spectra found.")
      
  datastart <- grep (sprintf ("^[[:blank:]]*##(%s)=", paste (.DATA.START, collapse = "|")), jdx) + 1
      # V 4.24 uses ##XYDATA=
      # V 5.00 uses ##DATA TABLE= ..., XYDATA
      # V 5.01 MPI Golm files use ##PEAK TABLE=

  if (length (datastart) == 0L) stop ("No data found: unsupported data type.")
      
  dataend <- grep ("^[[:blank:]]*##", jdx)
  dataend <- sapply (datastart, function (s) dataend [which (dataend > s)[1]]) - 1
  
  spcend <- grep ("^[[:blank:]]*##END=[[:blank:]]*$", jdx) - 1

  ## some checks
  stopifnot (length (datastart) >= length (hdrstart))
  stopifnot (length (datastart) == length (dataend))
  stopifnot (all (hdrstart < spcend))
  stopifnot (all (datastart < dataend))

  spc <- list ()

  for (s in seq_along (datastart)){
    ## look for header data
    hdr <- modifyList (header, .jdx.readhdr (jdx [hdrstart [s] : (datastart [s] - 1)]))

    if (! is.null (hdr$page) || ! is.null (hdr$ntuples))
        stop ("NTUPLES / PAGEs are not yet supported.")

    if (s == 1L) { ## file header may contain overall settings
        hdr <- modifyList (list (file = as.character (filename)), hdr)
        header <- hdr [! names (hdr) %in% .key2names (.DATA.START)]
      }
    
    ## evaluate data block

    if (grepl ("[A-DF-Za-df-z%@]", jdx[datastart [s]]))
        stop ("SQZ, DIF, and DIFDUP forms are not yet supported.")

    spc [[s]] <- switch (hdr$.format,
                         `(X++(Y..Y))`= .jdx.TABULAR.PAC  (hdr, jdx [datastart [s] : spcend [s]], ...),
                         `(XY..XY)`   = .jdx.TABULAR.AFFN (hdr, jdx [datastart [s] : spcend [s]], ...),

                         stop ("unknown JCAMP-DX data format: ", hdr$xydata)
                         )

    ## process according to header entries
    spc [[s]] <- .jdx.processhdr (spc [[s]], hdr, keys.hdr2data, ..., NA.symbols = NA.symbols)
  }

  if (length (spc) == 1L) 
    spc <- spc [[1]]
  else if (collapse.multi)
    spc <- collapse (spc)

  spc
}

### HEADER ------------------------------------------------------------------------------------------

.jdx.readhdr <- function (hdr){

  ## get rid of comments. JCAMP-DX comments start with $$ and go to the end of the line.
  hdr <- hdr [! grepl ("^[[:blank:]]*[$][$]", hdr)]
  hdr <- gsub ("([[:blank:]][$][$].*)$", "", hdr)

  ## now join lines that are not starting with ##KEY= with the KEYed line before
  nokey <- grep ("^[[:blank:]]*##.*=", hdr, invert = TRUE)
  if (length (nokey) > 0) {
    for (l in rev (nokey)) # these are few, so no optimization needed
        hdr [l - 1] <- paste (hdr [(l - 1) : l], collapse = " ")
    hdr <- hdr [-nokey]
  }

  names <-  .key2names (sub ("^[[:blank:]]*##(.*)=.*$", "\\1", hdr))
  
  hdr <- sub ("^[[:blank:]]*##.*=[[:blank:]]*(.*)[[:blank:]]*$", "\\1", hdr)
  hdr <-   gsub ("^[\"'[:blank:]]*([^\"'[:blank:]].*[^\"'[:blank:]])[\"'[:blank:]]*$", "\\1", hdr)
  i <- grepl ("^[[:blank:]]*[-]?[.[:digit:]]*[eE]?[-]?[.[:digit:]]*[[:blank:]]*$", hdr)
  hdr <- as.list (hdr)
  hdr [i] <- as.numeric (hdr [i])
  names (hdr) <- names

  ## e.g. Shimadzu does not always save XFACTOR and YFACTOR
  if (is.null (hdr$yfactor)) hdr$yfactor <- 1
  if (is.null (hdr$xfactor)) hdr$xfactor <- 1

  ## we treat XYDATA and PEAK TABLEs the same way
  format <- hdr [names (hdr) %in% .key2names (.DATA.START)]
  format <- format [! sapply (format, is.null)]
  if (length (format) != 1)
      stop ("contradicting format specification: please contact the maintainer (",
            maintainer ("hyperSpec"), 
            "supplying the file you just tried to load.")
  
  hdr$.format <- format [[1]]
      
  hdr
}

.jdx.processhdr <- function (spc, hdr, keys, ..., ytol = hdr$yfactor, NA.symbols){

  ## hdr$xfactor and $yfactor applied by individual reading functions

  ## check Y values
  miny <- min (spc@data$spc)
  if (! is.null (hdr$miny) && abs (hdr$miny - miny) > ytol)
      message (sprintf ("JDX file inconsistency: Minimum of spectrum != MINY: difference = %0.3g (%0.3g * YFACTOR)",
                        miny - hdr$miny,
                        (miny - hdr$miny) / hdr$yfactor))

  maxy <- max (spc@data$spc)
  if (! is.null (hdr$maxy) && abs (hdr$maxy - maxy) > ytol)
      message (sprintf ("JDX file inconsistency: Maximum of spectrum != MAXY: difference = %0.3g (%0.3g * YFACTOR)",
                        maxy - hdr$maxy,
                        (maxy - hdr$maxy) / hdr$yfactor))

  
  spc@label$.wavelength <- .jdx.xunits (hdr$xunits)
  spc@label$spc <- .jdx.yunits (hdr$yunits)

  ## CONCENTRATIONS
  if ("concentrations" %in% keys)
      spc <- .jdx.hdr.concentrations (spc, hdr, NA.symbols = NA.symbols)

  # delete header lines already processed
  hdr[c ("jcampdx", "xunits", "yunits", "xfactor", "yfactor", "firstx", "lastx", "npoints",
         "firsty", "xydata", "end", "deltax", "maxy", "miny", 
         "concentrations")] <- NULL
  if (is.character (keys))
      keys <- keys [keys %in% names (hdr)]
  hdr <- hdr [keys]

  if (length (hdr) > 0L)
      spc@data <- cbind (spc@data, hdr)
  
  spc
}

### DATA FORMATS ------------------------------------------------------------------------------------

.jdx.TABULAR.PAC <- function (hdr, data, ..., xtol = hdr$xfactor){

  ## regexp for numbers including scientific notation
  .PATTERN.number <- "[0-9]*[.]?[0-9]*([eE][+-]?[0-9]+)?"
  
  if (is.null (hdr$firstx))  stop ("##FIRSTX= missing.")
  if (is.null (hdr$lastx))   stop ("##LASTX= missing.")
  if (is.null (hdr$npoints)) stop ("##NPOINTS= missing.")
  
  wl <- seq (hdr$firstx, hdr$lastx, length.out = hdr$npoints)

  ## remove starting X
  y <- sub (paste0 ("^[[:blank:]]*", .PATTERN.number, "[[:blank:]]*(.*)$"), "\\2", data)

  ## add spaces between numbers if necessary
  y <- gsub ("([0-9.])([+-])", "\\1 \\2", y)

  y <- strsplit (y, "[[:blank:]]+")
  ny <- sapply (y, length)

  y <- as.numeric (unlist (y))


  if (length (y) != hdr$npoints)
      stop ("mismatch between ##NPOINTS and length of Y data.")

  ## X checkpoints
  x <- sub (paste0 ("^[[:blank:]]*(", .PATTERN.number, ")[[:blank:]]*.*$"), "\\1", data)
  x <- as.numeric (x) * hdr$xfactor
  diffx <- abs (wl [c (1, head (cumsum (ny) + 1, -1))] - x)
  if (any (diffx > xtol))
      message ("JDX file inconsistency: X axis differs from checkpoints. ",
               sprintf ("Maximum difference = %0.2g (%0.2g * XFACTOR)",
                        max (diffx), max (diffx) / hdr$xfactor))

  y <- y * hdr$yfactor

  new ("hyperSpec", spc = y, wavelength = wl)
}

.jdx.TABULAR.AFFN <- function (hdr, data, ...){
  #browser ()
  data <- strsplit (data, "[,;[:blank:]]+")
  data <- unlist (data)
  data <- matrix (as.numeric (data), nrow = 2)

  new ("hyperSpec", wavelength = data [1,] * hdr$xfactor, spc = data [2,]*hdr$yfactor)
}

### UNITS -------------------------------------------------------------------------------------------

.jdx.xunits <- function (xunits){
  if (is.null (xunits))
      NULL
  else 
      switch (tolower (xunits),
              `1/cm` = expression (tilde (nu) / cm^-1),
              micrometers = expression (`/` (lambda, micro * m)),
              nanometers  = expression (lambda / nm),
              seconds     = expression (t / s),
              xunits)
}

.jdx.yunits <- function (yunits){
  if (is.null (yunits))
      NULL

  else 
      switch (tolower (yunits),
              transmittance     = "T",
              reflectance       = "R",
              absorbance        = "A",
              `kubelka-munk`    = expression (`/` (1 - R^2, 2*R)),
              `arbitrary units` = "I / a.u.",
              yunits)
}

## HDR processing functions
.jdx.hdr.concentrations <- function (spc, hdr, NA.symbols){
  
  hdr <- strsplit (hdr$concentrations, "[)][[:blank:]]*[(]")[[1]]
  hdr [length (hdr)] <- gsub (")$", "", hdr [length (hdr)])
  if (hdr [1] == "(NCU")
      hdr <- hdr [-1]
  else
      message ("Unknown type of concentration specification in JDX file: ", hdr [1], ")")

  hdr <- simplify2array (strsplit (hdr, ","))
  hdr [hdr %in% NA.symbols] <- NA

  ## names
  N <- hdr [1,]
  N <- sub ("^([^[:alpha:]]*)", "", N)
  N <- sub ("([^[:alpha:]]*)$", "", N)
  N <- gsub ("([^[:alnum:]_-])", ".", N)
    
  ## concentrations
  C <- t (as.numeric (hdr [2,]))
  colnames (C) <- N
  C <- as.data.frame (C)
  spc@data <- cbind (spc@data, C)
    
  ## units
  U <- as.list (hdr [3,])
  names (U) <- N

  spc@label <- modifyList (spc@label, U)

  spc
}

## helpers
.DATA.START <- c ("XYDATA", "DATA TABLE", "PEAK TABLE")

.key2names <- function (key){
  gsub ("[[:blank:]_-]", "", tolower (key))
}
