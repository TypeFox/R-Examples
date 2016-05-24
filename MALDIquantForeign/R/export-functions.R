## Copyright 2013 Sebastian Gibb
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

.exportToFile <- function(x, file, type="auto", force=FALSE, ...) {
  if (file.exists(file) && !force) {
    stop("File already exists! Use ", sQuote("force=TRUE"), " to overwrite it.")
  }

  if (missing(type) || pmatch(tolower(type), "auto", nomatch=0,
                              duplicates.ok=FALSE)) {
    type <- .fileExtension(file)
  }

  i <- pmatch(tolower(type), exportFormats$type, nomatch=0, duplicates.ok=FALSE)

  if (i) {
    handler <- exportFormats$handler[i]
    return(do.call(handler, list(x=x, file=file, ...)))
  } else {
    stop("File type ", sQuote(type), " is not supported!")
  }
}

.exportToDir <- function(x, path, type, force=FALSE, ...) {
  if (!file.exists(path) && force) {
    dir.create(path, showWarnings=FALSE, recursive=TRUE)
  }

  if (!file.exists(path)) {
    stop("Directory ", sQuote(path), " doesn't exist!")
  }

  if (!file.info(path)$isdir) {
    stop(sQuote(path), " is no directory!")
  }

  ## stop if directory isn't writeable
  if (file.access(path, 2) != 0) {
    stop("No permissions to write into ", sQuote(path), "!")
  }

  i <- pmatch(tolower(type), exportFormats$type, nomatch=0, duplicates.ok=FALSE)

  if (i) {
    filenames <- .composeFilename(x, fileExtension=exportFormats$extension[i])
    filenames <- file.path(path, filenames)

    optArgs <- list(...)
    peaks <- list()

    if (hasArg(peaks)) {
      peaks <- optArgs$peaks
      optArgs$peaks <- NULL
    }

    for (i in seq(along=x)) {
      arguments <- list(x=x[[i]], file=filenames[i], type=type, force=force)
      arguments <- modifyList(arguments, optArgs)

      if (length(peaks)) {
        arguments$peaks <- peaks[[i]]
      }
      do.call(export, arguments)
    }
  } else {
    stop("File type ", sQuote(type), " is not supported!")
  }
  invisible()
}
