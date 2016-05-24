## Copyright 2012-2015 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of MALDIquantForeign for R and related languages.
##
## MALDIquantForeign is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## MALDIquantForeign is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with MALDIquantForeign. If not, see <http://www.gnu.org/licenses/>

#' Import files
#'
#' This function provides a general interface to import different file formats
#' into \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.
#'
#' @details
#' Specific import functions:
#' \tabular{ll}{
#'  txt \tab \code{\link[MALDIquantForeign]{importTxt}} \cr
#'  tab \tab \code{\link[MALDIquantForeign]{importTab}} \cr
#'  csv \tab \code{\link[MALDIquantForeign]{importCsv}} \cr
#'  fid \tab \code{\link[MALDIquantForeign]{importBrukerFlex}} \cr
#'  ciphergen \tab \code{\link[MALDIquantForeign]{importCiphergenXml}} \cr
#'  mzXML \tab \code{\link[MALDIquantForeign]{importMzXml}} \cr
#'  mzML \tab \code{\link[MALDIquantForeign]{importMzMl}} \cr
#'  imzML \tab \code{\link[MALDIquantForeign]{importImzMl}} \cr
#'  analyze \tab \code{\link[MALDIquantForeign]{importAnalyze}} \cr
#'  cdf \tab \code{\link[MALDIquantForeign]{importCdf}} \cr
#'  msd \tab \code{\link[MALDIquantForeign]{importMsd}} \cr
#' }
#'
#' \code{path}: In addition to the above mentioned file types the
#'  following (compressed) archives are supported, too:
#'  zip, tar, tar.gz, tar.bz2, tar.xz. The archives are uncompressed in a
#'  temporary directory. Afterwards the \code{\link[MALDIquantForeign]{import}}
#'  function is called (with \code{type="auto"}).
#'
#' \code{pattern}: Sometimes unusual file extensions are used (e.g.
#' \code{"*.xml"} for mzXML files). In this case a specific
#' \code{pattern} could be defined to import files with an unusual file
#' extension (e.g. \code{pattern="^.*\\.xml$"} to read all \code{*.xml}
#' files in a directory; see \code{\link[base]{regexp}} for details).
#'
#' \code{excludePattern}: Sometimes some files should be excluded. E.g.
#' to ignore additional aquired Bruker LIFT spectra
#' (MALDI-TOF/TOF; which are not support, yet) you could use
#' \code{excludePattern="([[:digit:]\\.]+)LIFT[\\\\/]1SRef"}.
#'
#' @param path \code{character}, path to directory or file which should be read
#'  in.
#' @param type \code{character}, file format. If \code{type} is set to
#'  \dQuote{auto} MALDIquant tries to detect the correct file type
#'  automatically. It often depends on the file extension
#'  (if \code{path} is a directory the most represented file extension is used;
#'  \code{pattern} argument is ignored).
#' @param pattern \code{character}, a regular expression to find files in a
#'  directory (see details).
#' @param excludePattern \code{character}, a regular expression to exclude
#'  files in a directory (see details).
#' @param removeEmptySpectra \code{logical}, should empty spectra excluded?
#' @param centroided \code{logical}, if \code{centroided=FALSE} (default)
#' the data are treated as not centroided and a list of
#' \code{\link[MALDIquant]{MassSpectrum-class}} objects is returned. Use
#' \code{centroided=TRUE} to assume centroided data and get a list of
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.
#' @param massRange \code{double}, limits of mass import (left/minimal mass,
#' right/maximal mass).
#' @param minIntensity \code{double}, minimal intensity to import.
#' @param mc.cores number of cores to use (default 1; only unix-based platforms
#' are supported, see
#' \code{\link[MALDIquantForeign]{MALDIquantForeign-parallel}} for details).
#' @param verbose \code{logical}, verbose output?
#' @param \ldots arguments to be passed to specific import functions.
#'
#' @return a \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects (depending on the
#' \code{centroided} argument).
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}}
#' \code{\link[MALDIquantForeign]{MALDIquantForeign-parallel}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}
#' @examples
#'
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' ## get example directory
#' exampleDirectory <- system.file("exampledata", package="MALDIquantForeign")
#'
#' ## import mzXML files
#' s <- import(exampleDirectory, type="mzXML")
#'
#' ## import tab delimited file with different file extension (default: *.tab)
#' s <- import(exampleDirectory, type="tab", pattern="^.*\\.txt")
#'
#' ## import single mzML file
#' s <- import(file.path(exampleDirectory, "tiny1.mzML1.1.mzML"))
#'
#' ## import gzipped csv file
#' s <- import(file.path(exampleDirectory, "compressed", "csv1.csv.gz"))
#'
#' @rdname import-functions
#' @export
import <- function(path, type="auto", pattern, excludePattern=NULL,
                   removeEmptySpectra=TRUE, centroided=FALSE, massRange=c(0, Inf),
                   minIntensity=0, mc.cores=1L, verbose=interactive(), ...) {

  ## download file if needed
  isUrl <- .isUrl(path)

  if (any(isUrl)) {
    path[isUrl] <- .download(path[isUrl], verbose=verbose)
    on.exit(.cleanupDownloadedTmpFiles())
  }

  ## file exists?
  isReadable <- file.exists(path) & file.access(path, mode=4) == 0

  if (any(!isReadable)) {
    stop(sQuote(path[!isReadable]), " doesn't exist or isn't readable!")
  }

  ## uncompress/unpack file if needed
  isCompressed <- .isPackedOrCompressed(path)

  if (any(isCompressed)) {
    path[isCompressed] <- .uncompress(path[isCompressed], verbose=verbose)
    on.exit(.cleanupUncompressedTmpFiles(), add=TRUE)
  }

  ## handle given file type
  i <- pmatch(tolower(type), c("auto", importFormats$type), nomatch=0,
              duplicates.ok=FALSE)-1

  if (i == -1) {
    stop("File type ", sQuote(type), " is not supported!")
  } else if (i == 0) {
    ## auto detect file type
    if (!missing(pattern)) {
      warning("User defined ", sQuote("pattern"), " is ignored in auto-mode.")
    }
    return(.importAuto(path=path, excludePattern=excludePattern,
                       removeEmptySpectra=removeEmptySpectra,
                       centroided=centroided, massRange=massRange,
                       minIntensity=minIntensity, verbose=verbose, ...))
  } else {
    ## user-defined file type
    if (missing(pattern)) {
      pattern <- importFormats$pattern[i]
    }
    handler <- get(importFormats$handler[i], mode="function")
    s <- unlist(MALDIquant:::.lapply(.files(path=path, pattern=pattern,
                                            excludePattern=excludePattern),
                                     handler, centroided=centroided,
                                     massRange=massRange,
                                     minIntensity=minIntensity,
                                     mc.cores=mc.cores,
                                     verbose=verbose, ...))
    if (is.null(s)) {
      stop("Import failed! Unsupported file type?")
    }

    if (removeEmptySpectra) {
      emptyIdx <- MALDIquant::findEmptyMassObjects(s)

      if (length(emptyIdx)) {
        .msg(verbose, "Remove ", length(emptyIdx), " empty spectra.")
        return(s[-emptyIdx])
      }
    }

    return(s)
  }
}

#' Import text files
#'
#' This function imports different text file formats
#' into \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.
#'
#' \code{importTab}, \code{importTxt} and \code{importCsv} use
#' \code{\link[utils]{read.table}} with different defaults.
#'
#' @param path \code{character}, path to directory or file which should be read
#'  in.
#' @param \ldots arguments to be passed to \code{\link[utils]{read.table}}.
#'
#' @return a \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects (depending on the
#' \code{centroided} argument).
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}},
#' \code{\link[utils]{read.table}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}
#' @examples
#'
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' ## get example directory
#' exampleDirectory <- system.file("exampledata", package="MALDIquantForeign")
#'
#' ## import txt files
#' s <- importTxt(exampleDirectory)
#'
#' ## import csv files
#' s <- importCsv(exampleDirectory)
#'
#' @rdname importTab-functions
#' @export
importTxt <- function(path, ...) {
  return(import(path=path, type="txt", ...))
}

#' @rdname importTab-functions
#' @export
importTab <- function(path, ...) {
  return(import(path=path, type="tab", ...))
}

#' @rdname importTab-functions
#' @export
importCsv <- function(path, ...) {
  return(import(path=path, type="csv", ...))
}

#' Import Bruker Daltonics *flex files
#'
#' This function imports files in Bruker Daltonics *flex-series file format
#' into \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.
#'
#' @param path \code{character}, path to directory or file which should be read
#'  in.
#' @param \ldots arguments to be passed to
#' \code{\link[readBrukerFlexData]{readBrukerFlexFile}}.
#'
#' @return a \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects (depending on the
#' \code{centroided} argument).
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}},
#' \code{\link[readBrukerFlexData]{readBrukerFlexFile}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}
#' @examples
#'
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' ## get example directory
#' exampleDirectory <- system.file("exampledata", package="MALDIquantForeign")
#'
#' s <- importBrukerFlex(exampleDirectory)
#'
#' @rdname importBrukerFlex-functions
#' @export
importBrukerFlex <- function(path, ...) {
  return(import(path=path, type="fid", ...))
}

#' Import mzXML files
#'
#' This function imports files in mzXML file format
#' into \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.
#'
#' @param path \code{character}, path to directory or file which should be read
#'  in.
#' @param \ldots arguments to be passed to
#' \code{\link[readMzXmlData]{readMzXmlFile}}.
#'
#' @return a \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects (depending on the
#' \code{centroided} argument).
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}},
#' \code{\link[readMzXmlData]{readMzXmlFile}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}, \cr
#' Definition of \code{mzXML} format:
#' \url{http://tools.proteomecenter.org/mzXMLschema.php}
#' @examples
#'
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' ## get example directory
#' exampleDirectory <- system.file("exampledata", package="MALDIquantForeign")
#'
#' ## import
#' s <- importMzXml(exampleDirectory)
#'
#' @rdname importMzXml-functions
#' @export
importMzXml <- function(path, ...) {
  return(import(path=path, type="mzxml", ...))
}

#' Import mzML files
#'
#' This function imports files in mzML file format
#' into \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.
#'
#' @param path \code{character}, path to directory or file which should be read
#'  in.
#' @param \ldots arguments to be passed to
#' \code{\link[MALDIquantForeign]{import}}.
#'
#' @return a \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects (depending on the
#' \code{centroided} argument).
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}, \cr
#' Definition of \code{mzML} format:
#' \url{http://www.psidev.info/mzml_1_0_0\%20}
#' @examples
#'
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' ## get example directory
#' exampleDirectory <- system.file("exampledata", package="MALDIquantForeign")
#'
#' ## import
#' s <- importMzMl(exampleDirectory)
#'
#' @rdname importMzMl-functions
#' @export
importMzMl <- function(path, ...) {
  return(import(path=path, type="mzml", ...))
}

#' Import imzML files
#'
#' This function imports files in imzML file format
#' into \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.
#'
#' @param path \code{character}, path to directory or file which should be read
#'  in.
#' @param coordinates \code{matrix}, 2 column matrix that contains the x- and
#'  y-coordinates for spectra that should be imported. Other spectra would be
#'  ignored.
#' @param \ldots arguments to be passed to
#' \code{\link[MALDIquantForeign]{import}}.
#'
#' @return a \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects (depending on the
#' \code{centroided} argument).
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}, \cr
#' Definition of \code{imzML} format:
#' \url{http://www.imzml.org/}
#' @examples
#'
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' ## get example directory
#' exampleDirectory <- system.file("exampledata", package="MALDIquantForeign")
#'
#' ## import
#' s <- importImzMl(file.path(exampleDirectory, "tiny_continuous.imzML"))
#'
#' ## import only spectra for pixel 1,1 and 2,1
#' s <- importImzMl(file.path(exampleDirectory, "tiny_continuous.imzML"),
#'                  coordinates = cbind(1:2, c(1, 1)))
#'
#' @rdname importImzMl-functions
#' @export
importImzMl <- function(path, coordinates=NULL, ...) {
  return(import(path=path, type="imzml", coordinates=coordinates, ...))
}

#' Import Ciphergen XML files
#'
#' This function imports files in Ciphergen XML file format
#' into \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.
#'
#' @param path \code{character}, path to directory or file which should be read
#'  in.
#' @param \ldots arguments to be passed to
#' \code{\link[MALDIquantForeign]{import}}.
#'
#' @return a \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects (depending on the
#' \code{centroided} argument).
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}
#' @examples
#'
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' ## get example directory
#' exampleDirectory <- system.file("exampledata", package="MALDIquantForeign")
#'
#' ## import
#' s <- importCiphergenXml(exampleDirectory)
#'
#' @rdname importCiphergenXml-functions
#' @export
importCiphergenXml <- function(path, ...) {
  return(import(path=path, type="ciphergen", ...))
}

#' Import Analyze 7.5 files
#'
#' This function imports files in Analyze 7.5 file format
#' into \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.
#'
#' @param path \code{character}, path to directory or file which should be read
#'  in.
#' @param \ldots arguments to be passed to
#' \code{\link[MALDIquantForeign]{import}}.
#'
#' @return a \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects (depending on the
#' \code{centroided} argument).
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/} \cr
#'  \url{http://www.grahamwideman.com/gw/brain/analyze/formatdoc.htm},
#'  \url{http://eeg.sourceforge.net/ANALYZE75.pdf}
#' @rdname importAnalyze-functions
#' @export
importAnalyze <- function(path, ...) {
  return(import(path=path, type="analyze", ...))
}

#' Import CDF files
#'
#' This function imports files in NetCDF file format
#' into \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.\cr
#' Please note that the \emph{RNetCDF} is needed.
#'
#' @param path \code{character}, path to directory or file which should be read
#'  in.
#' @param \ldots arguments to be passed to
#' \code{\link[MALDIquantForeign]{import}}.
#'
#' @return a \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects (depending on the
#' \code{centroided} argument).
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}
#' @examples
#'
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' ## get example directory
#' exampleDirectory <- system.file("exampledata", package="MALDIquantForeign")
#'
#' ## import
#' if (requireNamespace("RNetCDF", quietly=TRUE)) {
#'   s <- importCdf(exampleDirectory)
#' } else {
#'   message("You have to install the RNetCDF package to use importCdf.")
#' }
#'
#' @rdname importCdf-functions
#' @export
importCdf <- function(path, ...) {
  return(import(path=path, type="cdf", ...))
}

#' Import MSD files
#'
#' This function imports files in mMass MSD file format
#' into \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects.
#'
#' @param path \code{character}, path to directory or file which should be read
#'  in.
#' @param \ldots arguments to be passed to
#' \code{\link[MALDIquantForeign]{import}}.
#'
#' @return a \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} or
#' \code{\link[MALDIquant]{MassPeaks-class}} objects (depending on the
#' \code{centroided} argument).
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}, \cr
#' mMass homepage: \url{http://mmass.org/}
#' @examples
#'
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' ## get example directory
#' exampleDirectory <- system.file("exampledata", package="MALDIquantForeign")
#'
#' ## import
#' s <- importMsd(exampleDirectory)
#'
#' @rdname importMsd-functions
#' @export
importMsd <- function(path, ...) {
  return(import(path=path, type="msd", ...))
}
