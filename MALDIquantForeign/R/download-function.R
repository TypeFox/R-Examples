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

.isUrl <- function(x)grepl(pattern="^https?://|^ftp://", x=x)

.download <- function(url, destfile, verbose=FALSE, ...) {
  if (missing(destfile)) {
    pattern <- paste0(.withoutFileExtension(basename(url)), "_")
    fileext <- paste0(".", .fileExtension(url))
    tmpdir <- file.path(tempdir(), "MALDIquantForeign_download")

    if (!file.exists(tmpdir)) {
      dir.create(tmpdir, showWarnings=FALSE, recursive=TRUE)
    }

    destfile <- file.path(tmpdir, tempfile(pattern=pattern, tmpdir="",
                                           fileext=fileext))
  }

  .msg(verbose, "Downloading ", paste0(url, collapse=", ") , " to ",
       paste0(destfile, collapse=", "), ".")

  for (i in seq(along=url)) {
    if (download.file(url=url[i], destfile=destfile[i], quiet=!verbose,
                      mode="wb", ...)) {
      warning("Download of ", url[i], " failed!")
    }
  }

  destfile
}

.cleanupDownloadedTmpFiles <- function() {
  unlink(file.path(tempdir(), "MALDIquantForeign_download"), recursive=TRUE)
}
