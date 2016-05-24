## Copyright 2014 Sebastian Gibb
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

.testChecksum <- function(file, target, algo="sha1", ..., verbose=FALSE) {

  .msg(verbose, "Calculating ", algo, "-sum for ", sQuote(file), ": ",
       appendLF=FALSE)

  fileChecksum <- tolower(digest::digest(file, algo=algo, file=TRUE, ...))
  target <- tolower(target)

  .msg(verbose, fileChecksum)

  if (fileChecksum != target) {
    warning("Stored and calculated ", algo, " sums do not match ",
            "(stored: ", sQuote(target), ", calculated: ",
            sQuote(fileChecksum), ")!")
    return(FALSE)
  }

  return(TRUE)
}
