## Copyright 2012-2013 Sebastian Gibb
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

#' Export files
#'
#' This function provides a general interface to export
#' \code{\link[MALDIquant]{AbstractMassObject-class}} objects (e.g.
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}})
#' into different file formats.
#'
#' @details
#' Specific export functions:
#' \tabular{ll}{
#'  tab \tab \code{\link[MALDIquantForeign]{exportTab}} \cr
#'  csv \tab \code{\link[MALDIquantForeign]{exportCsv}} \cr
#'  imzML \tab \code{\link[MALDIquantForeign]{exportImzMl}} \cr
#'  msd \tab \code{\link[MALDIquantForeign]{exportMsd}} \cr
#'  mzML \tab \code{\link[MALDIquantForeign]{exportMzMl}} \cr
#' }
#'
#' @usage
#' \S4method{export}{AbstractMassObject}(x, file, type="auto", force=FALSE, \ldots)
#'
#' @param x a \code{\link[MALDIquant]{AbstractMassObject-class}} object or a
#'  \code{list} of \code{\link[MALDIquant]{AbstractMassObject-class}} objects.
#' @param file \code{character}, file name.
#' @param path \code{character}, path to directory in which the \code{list} of
#'  \code{\link[MALDIquant]{AbstractMassObject-class}} would be exported.
#' @param type \code{character}, file format. If \code{type} is set to
#'  \dQuote{auto} the file extension is used.
#' @param force \code{logical}, If \code{TRUE} the \code{file} would be
#'  overwritten or \code{path} would be created.
#' @param \ldots arguments to be passed to specific export functions.
#'
#' @seealso
#' \code{\link[MALDIquant]{MassPeaks-class}},
#' \code{\link[MALDIquant]{MassSpectrum-class}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}
#' @examples
#'
#' \dontrun{
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' s <- list(createMassSpectrum(mass=1:5, intensity=1:5),
#'           createMassSpectrum(mass=1:5, intensity=1:5))
#'
#' ## export a single spectrum
#' export(s[[1]], file="spectrum.csv")
#' ## identical to exportCsv(s[[1]], file="spectrum.csv")
#'
#' ## export a list of spectra
#' export(s, path="spectra", type="csv")
#' ## identical to exportCsv(s, path="spectra")
#' }
#'
#' @aliases export export,AbstractMassObject-method export,list-method
#' @rdname export-methods
#' @docType methods
#' @export
setMethod(f="export",
  signature=signature(x="AbstractMassObject"),
  definition=function(x, file, type="auto", force=FALSE, ...) {
  return(.exportToFile(x=x, file=file, type=type, force=force, ...))
})

#' @usage
#' \S4method{export}{list}(x, path, type, force=FALSE, \ldots)
#' @rdname export-methods
#' @export
setMethod(f="export",
  signature=signature(x="list"),
  definition=function(x, path, type, force=FALSE, ...) {

  stopifnot(MALDIquant:::.isMassObjectList(x))

  onefileSupport <- exportFormats$type[exportFormats$onefile]

  dots <- list(...)

  if (missing(path) && !is.null(dots$file)) {
    path <- dots$file
    dots$file <- NULL
  }

  isFile <- !isTRUE(file.info(path)$isdir) &&
            tolower(.fileExtension(path)) %in% onefileSupport

  if (isFile) {
    do.call(.exportToFile, modifyList(list(x=x, file=path, type=type, force=force),
                                      dots))
  } else {
    do.call(.exportToDir, modifyList(list(x=x, path=path, type=type, force=force),
                                     dots))
  }
})

#' Export to text files
#'
#' This function exports
#' \code{\link[MALDIquant]{AbstractMassObject-class}} objects (e.g.
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}})
#' into different text file formats.
#'
#' @details
#' \code{exportTab} and \code{exportCsv} use \code{\link[utils]{write.table}}
#' with different defaults (\code{sep="\t"} in \code{exportTab} and
#' \code{sep=","} in \code{exportCsv}).
#'
#' @usage
#' \S4method{exportTab}{AbstractMassObject}(x, file, force=FALSE, \ldots)
#'
#' @param x a \code{\link[MALDIquant]{AbstractMassObject-class}} object or a
#'  \code{list} of \code{\link[MALDIquant]{AbstractMassObject-class}} objects.
#' @param file \code{character}, file name.
#' @param path \code{character}, path to directory in which the \code{list} of
#'  \code{\link[MALDIquant]{AbstractMassObject-class}} would be exported.
#' @param force \code{logical}, If \code{TRUE} the \code{file} would be
#'  overwritten or \code{path} would be created.
#' @param \ldots arguments to be passed to \code{\link[utils]{write.table}}.
#'
#' @seealso
#' \code{\link[MALDIquant]{MassPeaks-class}},
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[utils]{write.table}}
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}
#' @examples
#'
#' \dontrun{
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' s <- list(createMassSpectrum(mass=1:5, intensity=1:5),
#'           createMassSpectrum(mass=1:5, intensity=1:5))
#'
#' ## export a single spectrum
#' exportTab(s[[1]], file="spectrum.tab")
#'
#' ## export a list of spectra and use ; as separator
#' exportCsv(s, path="spectra", sep=";", force=TRUE)
#' }
#'
#' @aliases exportTab exportTab,AbstractMassObject-method exportTab,list-method
#' exportCsv exportCsv,AbstractMassObject-method exportCsv,list-method
#' @rdname exportTab-methods
#' @docType methods
#' @export
setMethod(f="exportTab",
          signature=signature(x="AbstractMassObject"),
          definition=function(x, file, force=FALSE, ...) {
  export(x, file=file, type="tab", force=force, ...)
})

#' @usage
#' \S4method{exportTab}{list}(x, path, force=FALSE, \ldots)
#' @rdname exportTab-methods
#' @export
setMethod(f="exportTab",
          signature=signature(x="list"),
          definition=function(x, path, force=FALSE, ...) {
  export(x, path=path, type="tab", force=force, ...)
})

#' @usage
#' \S4method{exportCsv}{AbstractMassObject}(x, file, force=FALSE, \ldots)
#' @rdname exportTab-methods
#' @export
setMethod(f="exportCsv",
          signature=signature(x="AbstractMassObject"),
          definition=function(x, file, force=FALSE, ...) {
  export(x, file=file, type="csv", force=force, ...)
})

#' @usage
#' \S4method{exportCsv}{list}(x, path, force=FALSE, \ldots)
#' @rdname exportTab-methods
#' @export
setMethod(f="exportCsv",
          signature=signature(x="list"),
          definition=function(x, path, force=FALSE, ...) {
  export(x, path=path, type="csv", force=force, ...)
})

#' Export to MSD files
#'
#' This function exports
#' \code{\link[MALDIquant]{AbstractMassObject-class}} objects (e.g.
#' \code{\link[MALDIquant]{MassSpectrum-class}},
#' \code{\link[MALDIquant]{MassPeaks-class}})
#' into mMass MSD files.
#'
#' @usage
#' \S4method{exportMsd}{MassSpectrum}(x, file, force=FALSE, peaks, \ldots)
#'
#' @param x a \code{\link[MALDIquant]{MassSpectrum-class}} object or a
#'  \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} objects.
#' @param file \code{character}, file name.
#' @param path \code{character}, path to directory in which the \code{list} of
#'  \code{\link[MALDIquant]{AbstractMassObject-class}} would be exported.
#' @param peaks a \code{\link[MALDIquant]{MassPeaks-class}} object or a
#'  \code{list} of \code{\link[MALDIquant]{MassPeaks-class}} objects.
#' @param force \code{logical}, If \code{TRUE} the \code{file} would be
#'  overwritten or \code{path} would be created.
#' @param \ldots arguments to be passed to \code{\link[utils]{write.table}}.
#'
#' @seealso
#' \code{\link[MALDIquant]{MassPeaks-class}},
#' \code{\link[MALDIquant]{MassSpectrum-class}}
#'
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}, \cr
#' mMass homepage: \url{http://mmass.org/}
#' @examples
#'
#' \dontrun{
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' s <- list(createMassSpectrum(mass=1:5, intensity=1:5),
#'           createMassSpectrum(mass=1:5, intensity=1:5))
#' p <- list(createMassPeaks(mass=4:5, intensity=4:5, snr=1:2),
#'           createMassPeaks(mass=4:5, intensity=4:5, snr=1:2))
#'
#' ## export a single spectrum
#' exportMsd(s[[1]], file="spectrum.msd")
#'
#' ## export a single spectrum with corresponding peaks
#' exportMsd(s[[1]], file="spectrum.msd", peaks=p[[1]])
#'
#' ## export a list of spectra with corresponding peaks
#' exportMsd(s, path="spectra", peaks=p, force=TRUE)
#' }
#'
#' @aliases exportMsd exportMsd,MassSpectrum-method exportMsd,list-method
#' @rdname exportMsd-methods
#' @docType methods
#' @export
setMethod(f="exportMsd",
          signature=signature(x="MassSpectrum"),
          definition=function(x, file, force=FALSE, peaks, ...) {
  if (!missing(peaks)) {
    stopifnot(isMassPeaks(peaks))
    export(x, file=file, type="msd", force=force, peaks=peaks,  ...)
  } else {
    export(x, file=file, type="msd", force=force, ...)
  }
})

#' @usage
#' \S4method{exportMsd}{list}(x, path, force=FALSE, peaks, \ldots)
#' @rdname exportMsd-methods
#' @export
setMethod(f="exportMsd",
          signature=signature(x="list"),
          definition=function(x, path, force=FALSE, peaks, ...) {
  stopifnot(isMassSpectrumList(x))

  if (!missing(peaks)) {
    stopifnot(isMassPeaksList(peaks))
    export(x, path=path, type="msd", force=force, peaks=peaks,  ...)
  } else {
    export(x, path=path, type="msd", force=force, ...)
  }
})

#' Export to mzML files
#'
#' This function exports
#' \code{\link[MALDIquant]{MassSpectrum-class}} objects into mzML files.
#'
#' @usage
#' \S4method{exportMzMl}{MassSpectrum}(x, file, force=FALSE, \ldots)
#'
#' @param x a \code{\link[MALDIquant]{MassSpectrum-class}} object or a
#'  \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} objects.
#' @param file \code{character}, file name.
#' @param path \code{character}, path to directory in which the \code{list} of
#'  \code{\link[MALDIquant]{MassSpectrum-class}} would be exported. If
#'  \code{path} is a single filename all spectra will be exported to a single
#'  mzML file.
#' @param force \code{logical}, If \code{TRUE} the \code{file} would be
#'  overwritten or \code{path} would be created.
#' @param \ldots arguments to be passed to internal functions.
#'
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}}
#'
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}, \cr
#' HUPO Proteomics Standards Inititative mzML 1.1.0 Specification:
#' \url{http://www.psidev.info/mzml_1_0_0}
#' @examples
#'
#' \dontrun{
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' s <- list(createMassSpectrum(mass=1:5, intensity=1:5),
#'           createMassSpectrum(mass=1:5, intensity=1:5))
#'
#' ## export a single spectrum
#' exportMzMl(s[[1]], file="spectrum.mzML")
#'
#' ## export a list of spectra
#' exportMzMl(s, path="spectra.mzML")
#' }
#'
#' @aliases exportMzMl exportMzMl,MassSpectrum-method exportMzMl,list-method
#' @rdname exportMzMl-methods
#' @docType methods
#' @export
setMethod(f="exportMzMl",
          signature=signature(x="MassSpectrum"),
          definition=function(x, file, force=FALSE, ...) {
  export(x, file=file, type="mzml", force=force, ...)
})

#' @usage
#' \S4method{exportMzMl}{list}(x, path, force=FALSE, \ldots)
#' @rdname exportMzMl-methods
#' @export
setMethod(f="exportMzMl",
          signature=signature(x="list"),
          definition=function(x, path, force=FALSE, ...) {
  export(x, path=path, type="mzml", force=force, ...)
})

#' Export to imzML files
#'
#' This function exports
#' \code{\link[MALDIquant]{MassSpectrum-class}} objects into imzML files.
#'
#' @usage
#' \S4method{exportImzMl}{MassSpectrum}(x, file, force=FALSE, processed=TRUE,
#' coordinates=NULL, pixelSize=c(100, 100), \ldots)
#'
#' @param x a \code{\link[MALDIquant]{MassSpectrum-class}} object or a
#'  \code{list} of \code{\link[MALDIquant]{MassSpectrum-class}} objects.
#' @param file \code{character}, file name.
#' @param path \code{character}, path to directory in which the \code{list} of
#'  \code{\link[MALDIquant]{MassSpectrum-class}} would be exported. If
#'  \code{path} is a single filename all spectra will be exported to a single
#'  imzML file.
#' @param force \code{logical}, If \code{TRUE} the \code{file} would be
#'  overwritten or \code{path} would be created.
#' @param processed \code{logical}, If \code{TRUE} (default) the spectra will
#' be saved in processed mode (means mass and intensity is stored for each
#' spectra separately in contrast to continuous mode where the mass is stored
#' only for one spectrum).
#' @param coordinates \code{matrix}, 2 column matrix that contains the x- and
#'  y-coordinates for the spectra.
#' @param pixelSize \code{numeric}, a vector of length 2 that contains the x and
#' y pixel size in micrometers (default: \code{c(100, 100)}).
#' @param \ldots arguments to be passed to internal functions.
#'
#' @seealso
#' \code{\link[MALDIquant]{MassSpectrum-class}}
#'
#' @author Sebastian Gibb
#' @references \url{http://strimmerlab.org/software/maldiquant/}
#'
#' Schramm T, Hester A, Klinkert I, Both J-P, Heeren RMA, Brunelle A,
#' Laprevote O, Desbenoit N, Robbe M-F, Stoeckli M, Spengler B, Roempp A
#' (2012)\cr
#' imzML - A common data format for the flexible exchange and processing of mass
#' spectrometry imaging data.\cr
#' Journal of Proteomics 75 (16):5106-5110. \cr
#' \url{http://dx.doi.org/10.1016/j.jprot.2012.07.026}
#' @examples
#'
#' \dontrun{
#' library("MALDIquant")
#' library("MALDIquantForeign")
#'
#' s <- list(createMassSpectrum(mass=1:5, intensity=1:5),
#'           createMassSpectrum(mass=1:5, intensity=1:5))
#'
#' ## export a list of spectra
#' exportImzMl(s, path="processed.imzML", coordinates=cbind(x=1:2, y=c(1, 1)))
#' }
#'
#' @aliases exportImzMl exportImzMl,MassSpectrum-method exportImzMl,list-method
#' @rdname exportImzMl-methods
#' @docType methods
#' @export
setMethod(f="exportImzMl",
          signature=signature(x="MassSpectrum"),
          definition=function(x, file, force=FALSE,
                              processed=TRUE, coordinates=NULL,
                              pixelSize=c(100, 100), ...) {
  export(x, file=file, type="imzml", force=force,
         processed=processed, coordinates=coordinates, pixelSize=pixelSize, ...)
})

#' @usage
#' \S4method{exportImzMl}{list}(x, path, force=FALSE, processed=TRUE,
#' coordinates=NULL, pixelSize=c(100, 100), \ldots)
#' @rdname exportImzMl-methods
#' @export
setMethod(f="exportImzMl",
          signature=signature(x="list"),
          definition=function(x, path, force=FALSE,
                              processed=TRUE, coordinates=NULL,
                              pixelSize=c(100, 100), ...) {
  export(x, path=path, type="imzml", force=force,
         processed=processed, coordinates=coordinates, pixelSize=pixelSize, ...)
})
