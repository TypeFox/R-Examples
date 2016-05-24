## Copyright 2010-2014 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of readBrukerFlexData for R and related languages.
##
## readBrukerFlexData is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## readBrukerFlexData is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with readBrukerFlexData. If not, see <http://www.gnu.org/licenses/>

#' Creates sample name.
#'
#' This function guesses the name of current spot from filename.
#'
#' @note WARNING: if the 4th upper directory hasn't an unique name you will get
#'  equal names for your spot list
#'
#'  e.g. the following will create a list with 4 elements but
#'  only 2 unique spot names (2-100kDa:0_A1 and 2-100kDa:0_B1) \cr
#' ./Run1/2-100kDa/0_A1/1/1SLin/fid \cr
#' ./Run1/2-100kDa/0_B1/1/1SLin/fid \cr
#' ./Run2/2-100kDa/0_A1/1/1SLin/fid \cr
#' ./Run2/2-100kDa/0_B1/1/1SLin/fid \cr
#'
#' @param fidFile \code{character}, path to fid file
#'  (e.g. \sQuote{Pankreas_HB_L_061019_A10/0_a19/1/1SLin/fid})
#'
#' @return
#'  \code{character}, sample name: e.g \sQuote{Pankreas_HB_L_061019_A10}
#'
#' @seealso
#'  \code{\link[readBrukerFlexData]{.readAcquFile}},
#'  \code{\link[readBrukerFlexData]{readBrukerFlexFile}}
#' @keywords internal
#' @rdname sampleName
#'
.sampleName <- function(fidFile) {
  # regular expression for directory separator (on unix: /+, on windows \+)
  # sadly .Platform$file.sep == "/" on both
  fidFile <- chartr(old="\\", new="/", x=fidFile)

  # create array of directories (each element == one directory)
  dirs <- strsplit(x=fidFile, split="/")[[1L]]

  numDirs <- length(dirs)

  sampleName <- NA

  ## old FlexAnalysis seems to have the following directories
  ## 0_L20_1SLin/fid
  ## vs the more recent FlexAnalysis versions use
  ## 0/L20/1SLin/fid
  isShortPath <- isTRUE(numDirs > 2 &&
                        grepl("[0-9]+_[A-z][0-9]+_[0-9][A-z]+$",
                              dirs[numDirs - 1L]))
  if (isShortPath) {
    sampleName <- dirs[numDirs-2L]
  } else if (numDirs > 4L ) {
    sampleName <- dirs[numDirs-4L]

    # -, : or something like that causes errors in names()
    # TODO: use make.names in future releases?
    sampleName <- gsub(pattern="[[:punct:]]|[[:space:]]", replacement="_",
                       x=sampleName)
  }

  return(sampleName)
}

