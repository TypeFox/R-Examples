####################################################################################################
###
###  read.ENVI - read ENVI files, missing header files may be replaced by list in parameter header 
###
###  * read.ENVI.Nicolet for ENVI files written by Nicolet spectrometers 
###  * adapted from caTools read.ENVI
###
###  Time-stamp: <Claudia Beleites on Saturday, 2011-02-05 at 19:19:00 on cb>
###
####################################################################################################

### some general helper functions ..................................................................

###-----------------------------------------------------------------------------
###
### split.line - split line into list of key-value pairs
###
###

split.line <- function (x, separator, trim.blank = TRUE) {
  tmp <- regexpr (separator, x)

  key   <- substr (x, 1, tmp - 1)
  value <- substr (x, tmp + 1, nchar (x))

  if (trim.blank){
    blank.pattern <- "^[[:blank:]]*([^[:blank:]]+.*[^[:blank:]]+)[[:blank:]]*$"
    key <- sub (blank.pattern, "\\1", key)
    value <- sub (blank.pattern, "\\1", value)
  }

  value <- as.list (value)
  names (value) <- key

  value
}


### some ENVI-specific helper functions .............................................................


.find.ENVI.header  <- function (file, headerfilename) {
  if (is.null (headerfilename)) {
    headerfilename <- paste (dirname (file), sub ("[.][^.]+$", ".*", basename (file)), sep = "/")
    tmp <- Sys.glob (headerfilename)
    headerfilename <- tmp [! grepl (file, tmp)]

    if (length (headerfilename) > 1) {
      headerfilename <- headerfilename [grepl ("[.][hH][dD][rR]$", headerfilename)]
      if (length (headerfilename == 1))
        message (".read.ENVI.header: Guessing header file name ", headerfilename)
    }
    
    if (length (headerfilename) != 1)
      stop ("Cannot guess header file name")
    else
      message (".read.ENVI.header: Guessing header file name ", headerfilename)
  }
    
  if (!file.exists(headerfilename)) 
    stop("Could not open header file: ", headerfilename)

  headerfilename
}

# ...................................................................................................

.read.ENVI.split.header <- function (header) {
  ## check ENVI at beginning of file
  if (!grepl("ENVI", header[1])) 
    stop("Not an ENVI header (ENVI keyword missing)")
  else
    header <- header [-1]

  ## remove curly braces and put multi-line key-value-pairs into one line
  header <- gsub("\\{([^}]*)\\}", "\\1", header)

  l <- grep("\\{", header)
  r <- grep("\\}", header)
  
  if (length(l) != length(r)) 
    stop("Error matching curly braces in header (differing numbers).")

  if (any(r <= l)) 
    stop("Mismatch of curly braces in header.")

  header[l] <- sub("\\{", "", header[l])
  header[r] <- sub("\\}", "", header[r])

  for (i in rev(seq_along(l))) 
    header <- c(header[seq_len(l[i] - 1)],
                paste(header[l[i]:r[i]], collapse = " "),
                header[-seq_len(r[i])])

  ## split key = value constructs into list with keys as names
  header <- sapply (header, split.line, "=", USE.NAMES = FALSE)
  names (header) <- tolower (names (header))

  ## some values are numeric
  tmp <- names (header) %in% c("samples", "lines", "bands", "data type", "header offset")
  header [tmp] <- lapply (header [tmp], as.numeric)

  header
}

### .................................................................................................

.read.ENVI.bin <- function (file, header) {
  if (any (is.null (header [c("samples", "lines", "bands", "data type")]) ||
           is.na   (header [c("samples", "lines", "bands", "data type")]) ))
    stop("Error in ENVI header (required entry missing or incorrect)\n header: ",
         paste (names (header), " = ", header, collapse = ", "))

  if (header$samples <= 0)
    stop("Error in ENVI header: incorrect data size (", header$samples, ")")
  if (header$lines <= 0)
    stop("Error in ENVI header: incorrect data size (", header$lines, ")")
  if (header$bands <= 0)
    stop("Error in ENVI header: incorrect data size (", header$bands, ")")
  
  if (!(header$`data type` %in% c(1 : 5, 9, 12))) 
    stop("Error in ENVI header: data type incorrect or unsupported (", header$`data type`,")")

  if (is.null (header$`byte order`)){
    header$`byte order` <- .Platform$endian
    message (".read.ENVI.bin: 'byte order' not given => Guessing '",
             .Platform$endian, "'\n", sep = '')
  }
  if (! header$`byte order` %in% c ("big", "little", "swap")) {
    header$`byte order` <- as.numeric (header$`byte order`)
    if (! header$`byte order` %in% 0 : 1) {
      header$`byte order` <- .Platform$endian
      warning ("byte order incorrect. Guessing '", .Platform$endian, "'")
    } else if (header$`byte order` == 0)
      header$`byte order` <- "little"
    else 
      header$`byte order` <- "big"
  }

  n <- header$samples * header$lines * header$bands

  if (!file.exists(file)) 
    stop("Could not open binary file: ", file)

  f <- file (file, "rb")
  if (! is.null (header$`header offset`)) 
    readBin(f, raw(), n = header$`header offset`)
  
  switch(header$`data type`,
         spc <- readBin(f, integer(), n = n, size =  1, signed = FALSE),
         spc <- readBin(f, integer(), n = n, size =  2, endian = header$`byte order`),
         spc <- readBin(f, integer(), n = n, size =  4, endian = header$`byte order`),
         spc <- readBin(f, double(),  n = n, size =  4, endian = header$`byte order`),
         spc <- readBin(f, double(),  n = n, size =  8, endian = header$`byte order`),
         , # 6 unused
         , # 7 unused
         , # 8 unused
         spc <- readBin(f, complex(), n = n, size = 16, endian = header$`byte order`),
         , # 10 unused
         , # 11 unused
         spc <- readBin(f, integer(), n = n, size =  2, endian = header$`byte order`, signed = FALSE)
         )
  
  close(f)

  if (is.null (header$interleave))
    header$interleave <- "bsq"    # de
  
  switch (tolower (header$interleave),
          bil = {
            dim (spc) <- c(header$samples, header$bands, header$lines);
            spc <- aperm(spc, c(3, 1, 2))
          },
          bip = {
            dim (spc) <- c(header$bands, header$samples, header$lines);
            spc <- aperm(spc, c(3, 2, 1))
          },
          bsq = {
            dim (spc) <- c(header$samples, header$lines, header$bands);
            spc <- aperm(spc, c(2, 1, 3))
          },
          stop ("Unknown interleave (",
                header$interleave,
                ", should be one of 'bsq', 'bil', 'bip')")
          )

  dim (spc) <- c (header$samples * header$lines, header$bands)

  spc
}

# ..................................................................................................



##' Import of ENVI data as hyperSpec object
##' This function allows ENVI data import as \code{hyperSpec} object.
##' 
##' \code{read.ENVI.Nicolet} should be a good starting point for writing custom
##' wrappers for \code{read.ENVI} that take into account your manufacturer's
##' special entries in the header file.
##' 
##' ENVI data usually consists of two files, an ASCII header and a binary data
##' file. The header contains all information necessary for correctly reading
##' the binary file.
##' 
##' I experienced missing header files (or rather: header files without any
##' contents) produced by Bruker Opus' ENVI export.
##' 
##' In this case the necessary information can be given as a list in parameter
##' \code{header} instead. The elements of header are then:
##' 
##' \tabular{lll}{
##' \code{header\$}         \tab values        \tab meaning\cr
##' \code{samples}          \tab integer       \tab no of columns / spectra in x direction\cr
##' \code{lines}            \tab integer       \tab no of lines / spectra in y direction\cr
##' \code{bands}            \tab integer       \tab no of wavelengths / data points per spectrum\cr
##' \code{`data type`}      \tab               \tab format of the binary file\cr
##'                         \tab 1             \tab 1 byte unsigned integer \cr
##'                         \tab 2             \tab 2 byte signed integer \cr
##'                         \tab 3             \tab 4 byte signed integer \cr
##'                         \tab 4             \tab 4 byte float \cr
##'                         \tab 5             \tab 8 byte double \cr
##'                         \tab 9             \tab 16 (2 x 8) byte complex double \cr
##'                         \tab 12            \tab 2 byte unsigned integer \cr
##'  \code{`header offset`} \tab integer       \tab number of bytes to skip before binary data starts\cr
##'  \code{interleave}      \tab               \tab directions of the data cube \cr
##'                         \tab "BSQ"         \tab band sequential (indexing: [sample, line, band])\cr
##'                         \tab "BIL"         \tab band interleave by line (indexing: [sample, line, band])\cr
##'                         \tab "BIP"         \tab band interleave by pixel (indexing: [band, line, sample])\cr
##'  \code{`byte order`}    \tab 0 or "little" \tab little endian \cr
##'                         \tab 1 or "big"    \tab big endian \cr
##'                         \tab "swap"        \tab swap byte order
##' }
##' 
##' Some more information that is not provided by the ENVI files may be given:
##' 
##' Wavelength axis and axis labels in the respective parameters. For more
##' information, see \code{\link[hyperSpec]{initialize}}.
##' 
##' The spatial information is by default a sequence from 0 to
##' \code{header$samples - 1} and \code{header$lines - 1}, respectively.
##' \code{x} and \code{y} give offset of the first spectrum and step size.
##' 
##' Thus, the object's \code{$x} colum is: \code{(0 : header$samples - 1) * x
##' [2] + x [1]}.  The \code{$y} colum is calculated analogously.
##' @aliases read.ENVI read.ENVI.Nicolet
##' @param file complete name of the binary file
##' @param headerfile name of the ASCII header file. If \code{NULL}, the name
##'   of the header file is guessed by looking for a second file with the same
##'   basename but different suffix as \code{file}.
##' @param header list with the respective information, see details.
##' @param x,y vectors of form c(offset, step size) for the position vectors,
##'   see details.
##' @param wavelength,label lists that overwrite the respective information
##'   from the ENVI header file. These data is then handed to
##'   \code{\link[hyperSpec]{initialize}}
##' @param keys.hdr2data determines which fields of the header file should be put into the extra
##' data. Defaults to none.
##' 
##' To specify certain entries, give character vectors containing the lowercase
##'   names of the header file entries.
##' @return a \code{hyperSpec} object
##' @author C. Beleites, testing for the Nicolet files C. Dicko
##' @seealso \code{\link[caTools]{read.ENVI}}
##' 
##' \code{\link[hyperSpec]{textio}}
##' @references This function was adapted from
##'   \code{\link[caTools]{read.ENVI}}:
##' 
##' Jarek Tuszynski (2008). caTools: Tools: moving window statistics, GIF,
##'   Base64, ROC AUC, etc.. R package version 1.9.
##' @rdname readENVI
##' @export
##' @keywords IO file
read.ENVI <- function (file = stop ("read.ENVI: file name needed"), headerfile = NULL, 
							  header = list (), 
							  keys.hdr2data = FALSE, 
							  x = 0 : 1, y = x, 
							  wavelength = NULL, label = list ()) {
  force (y)

  if (! file.exists (file))
	  stop ("File not found:", file)

  if (! is.list (header)) # catch a common pitfall
    if (is.character (header))
      stop ("header must be a list of parameters. Did you mean headerfile instead?")
    else
      stop ("header must be a list of parameters.")
				
  if (is.null (headerfile))
  	headerfile <- .find.ENVI.header (file, headerfile)
  
  tmp <- readLines (headerfile)
  tmp <- .read.ENVI.split.header (tmp)
  header <- modifyList (tmp, header)  

  ## _no_ capital letters here: .read.ENVI.split.header produces lowercase names
  recognized.keywords <- c("samples", "lines", "bands", "data type", "header offset", 
                           "interleave", "byte order", "wavelength")

  ## read the binary file
  spc <- .read.ENVI.bin (file, header)

  ## wavelength should contain the mean wavelength of the respective band
  if (! is.null (header$wavelength)) {
    header$wavelength <- as.numeric (unlist (strsplit (header$wavelength, "[,;[:blank:]]+")))

    if (! any (is.na (header$wavelength)) && is.null (wavelength))
      wavelength <- header$wavelength
  } 
  
  ## set up spatial coordinates
  x <- rep (seq (0, header$samples - 1), each = header$lines)   * x [2] + x [1]
  y <- rep (seq (0, header$lines   - 1),        header$samples) * y [2] + y [1]
  
  ## header lines => extra data columns 
  extra.data <- header [keys.hdr2data]

  if (.options$gc) gc ()
  
  if (length (extra.data) > 0) {
	  extra.data <- lapply (extra.data, rep, length.out = length (x))
	  data <- data.frame (x = x, y = y, extra.data)
  } else {
	  data <- data.frame (x = x, y = y)
  }
  
  if (.options$gc) gc ()

  ## finally put together the hyperSpec object
  spc <- new ("hyperSpec", data = data, spc = spc, wavelength = wavelength, labels = label)
  
  ## consistent file import behaviour across import functions
  .fileio.optional (spc, file)
}

