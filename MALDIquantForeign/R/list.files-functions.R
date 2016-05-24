## Copyright 2012 Sebastian Gibb
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

.list.files <- function(path, pattern, excludePattern=NULL, recursive=TRUE,
                        ignore.case=TRUE) {
  files <- list.files(path=path, pattern=pattern, recursive=recursive,
                      ignore.case=ignore.case, full.names=TRUE)

  if (!is.null(excludePattern)) {
    isExcluded <- grepl(pattern=excludePattern, x=files,
                        ignore.case=ignore.case)
    files <- files[!isExcluded]
  }

  normalizePath(files)
}

.files <- function(path, pattern, excludePattern=NULL, ignore.case=TRUE, ...) {
  isDir <- file.info(path)$isdir

  files <- normalizePath(path[!isDir])

  isMatching <- unlist(regexpr(pattern=pattern, text=basename(files),
                               ignore.case=ignore.case)) != -1

  files <- files[isMatching]
  files <- c(files, .list.files(path=path[isDir], pattern=pattern,
                                excludePattern=excludePattern,
                                ignore.case=ignore.case, ...))
  unique(files)
}
