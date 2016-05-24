## Copyright 2012-2014 Sebastian Gibb
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

.importTab <- function(file, centroided=FALSE, massRange=c(0L, Inf),
                       minIntensity=0L, skip=0L,
                       sep=.autoSep(file, skip=skip),
                       header=.autoHeader(file, sep=sep, skip=skip),
                       verbose=FALSE, ...) {
  ## load ms file
  s <- read.table(file=file, header=header, sep=sep, skip=skip,
                  stringsAsFactors=FALSE, ...)

  list(.createMassObject(mass=s[, 1L], intensity=s[, 2L],
                         metaData=list(file=file),
                         centroided=centroided,
                         massRange=massRange,
                         minIntensity=minIntensity,
                         verbose=verbose))
}

.importCsv <- function(file, centroided=FALSE, massRange=c(0L, Inf),
                       minIntensity=0L, skip=0L,
                       sep=.autoSep(file, skip=skip),
                       header=.autoHeader(file, sep=sep, skip=skip),
                       verbose=FALSE, ...) {

  .importTab(file=file, centroided=centroided, massRange=massRange,
             minIntensity=minIntensity, skip=skip, sep=sep,
             header=header, verbose=verbose, ...)
}

.autoHeader <- function(file, sep="\t", skip=0L) {
  l <- tail(readLines(file, n=skip+1L), 1L)
  l <- gsub(pattern='[\\\\"]*', replacement="", x=l)
  l <- strsplit(l, split=sep)[[1L]][1L]
  !is.numeric(type.convert(l, as.is=TRUE))
}

.autoSep <- function(file, sep=c(",", ";", "\t", " "), skip=0L) {
  l <- tail(readLines(file, n=skip+1L), 1L)
  pattern <- paste0(".+", sep, ".+")
  i <- vapply(pattern, function(x) {
    g <- gregexpr(pattern=x, text=l)[[1L]]
    return(all(g > 0L) & length(g) == 1L)
  }, logical(1L))

  if (any(i)) {
    return(sep[which(i)[1L]])  ## return only first match
  } else {
    return(sep[1L])
  }
}
