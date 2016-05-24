##    This file is part of sciViews.
##
##    sciViews is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    sciViews is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with sciViews.  If not, see <http://www.gnu.org/licenses/>.
##

makeTestListFromExamples <- function(packageName, manFilesDir, skipFailing=FALSE) {
  manPageFiles <- list.files(manFilesDir, pattern="\\.Rd$")
  manPages <- sapply(manPageFiles, function(filename) {
    lines <- readLines(paste(manFilesDir, filename, sep="/"))
    lines <- lines[grep("^\\\\name[ ]*\\{(.*)\\}", lines)]
    sub("^\\\\name[ ]*\\{(.*)\\}", lines, replacement="\\1")
  })
  manPages <- manPages[manPages != paste(packageName, "package", sep="-")]
  names(manPages) <- manPages

  lapply(manPages, function(x) {
    testCall <- call("example", x, packageName)
    if (skipFailing)
      onFailure <- function(w) { DEACTIVATED(paste(w, collapse="\n")); checkTrue(FALSE); }
    else
      onFailure <- function(w) checkIdentical(NULL, w)
    result <- svTest(function() {
                     tryCatch(withCallingHandlers({ eval(testCall); checkTrue(TRUE); },
                                                  warning=function(w) onFailure(w)),
                              error=function(w) onFailure(w))
                   })
    attr(result, 'unit') <- 'rcheck.examples'
    result
  })
}
