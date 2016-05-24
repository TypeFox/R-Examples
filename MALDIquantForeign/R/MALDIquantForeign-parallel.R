## Copyright 2015 Sebastian Gibb
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

#' Parallel Support in Package \pkg{MALDIquantForeign}
#'
#' \code{\link[MALDIquantForeign]{MALDIquantForeign-package}} offers multi-core
#' support using \code{\link[parallel]{mclapply}} and
#' \code{\link[parallel]{mcmapply}}. This approach is limited to unix-based
#' platforms.
#'
#' Please note that not all import functions benfit from parallelisation. The
#' current implementation is limited to run the parallelisation over different
#' files. That's why only imports of multiple files could be run on multiple
#' cores. E.g. a single mzML file containing 4 spectra would always be
#' read on a single core. In contrast 4 mzML files each containing just one
#' spectra could be read in using 4 cores.
#'
#' The improvement in the runtime depends on the amount of data to read, the
#' proportion of parsing/decoding of the data, the amount of memory and the
#' speed of the hard disk.
#'
#' Please note: It is possible that using parallelisation results in a worse
#' runtime!
#'
#' @name MALDIquantForeign-parallel
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @references \url{http://strimmerlab.org/software/maldiquant/}
#' @keywords misc
#' @seealso
#' \code{\link[MALDIquant]{MALDIquant-parallel}},
#' \code{\link[parallel]{mclapply}},
#' \code{\link[parallel]{mcmapply}}
#'
#' @examples
#'  ## load packages
#'  library("MALDIquant")
#'  library("MALDIquantForeign")
#'
#'  exampleDirectory <- system.file("exampledata", package="MALDIquantForeign")
#'
#'  ## run single-core import
#'  print(system.time(
#'    s1 <- importMzMl(exampleDirectory, centroided=TRUE, verbose=FALSE)
#'  ))
#'
#'  if(.Platform$OS.type == "unix") {
#'    ## run multi-core import
#'    ## (because the example spectra are very small (just 5 data points) the
#'    ## multi-core solution is slower on most systems)
#'    print(system.time(
#'      s2 <- importMzMl(exampleDirectory, centroided=TRUE, mc.cores=2,
#'                       verbose=FALSE)
#'    ))
#'    stopifnot(all.equal(s1, s2))
#'  }
#' @rdname MALDIquantForeign-parallel

NULL
