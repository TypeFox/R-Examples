## rmatio, a R interface to the C library matio, MAT File I/O Library.
## Copyright (C) 2013-2014  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## rmatio is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' Writes the values in a list to a mat-file.
##'
##' Writes the values in the list to a mat-file. All values in the
##' list must have unique names.
##' @note
##' \itemize{
##'   \item A vector is saved as a \code{1 x length} array
##'
##'   \item Support for writing a sparse matrix of type 'dgCMatrix' or 'lgCMatrix'
##'     to file
##' }
##' @rdname write.mat-methods
##' @aliases write.mat
##' @aliases write.mat-methods
##' @aliases write.mat,list-method
##' @docType methods
##' @title Write Matlab file
##' @param object The \code{object} to write.
##' @param filename The MAT file to write.
##' @param compression Use compression when writing variables. Defaults to TRUE.
##' @param version MAT file version to create. Currently only support
##' for Matlab level-5 file (MAT5) from rmatio package.
##' @return invisible NULL
##' @keywords methods
##' @include have_lib.r
##' @export
##' @useDynLib rmatio
##' @author Stefan Widgren
##' @examples
##' \dontrun{
##' filename <- tempfile(fileext = ".mat")
##'
##' ## Example how to read and write an integer vector with rmatio
##' write.mat(list(a=1:5), filename=filename)
##' a <- as.integer(read.mat(filename)[["a"]])
##'
##' stopifnot(identical(a, 1:5))
##'
##' unlink(filename)
##'
##' ## Read a compressed version 5 MAT file
##' m <- read.mat(system.file("extdata/matio_test_cases_compressed_le.mat",
##'                           package='rmatio'))
##'
##' ## Write an uncompressed version 5 MAT file
##' write.mat(m, filename="test-uncompressed.mat", compression=FALSE, version="MAT5")
##'
##' ## Write a compressed version 5 MAT file
##' write.mat(m, filename="test-compressed.mat", compression=TRUE, version="MAT5")
##'
##' ## Check that the content of the files are identical
##' identical(read.mat("test-uncompressed.mat"), read.mat("test-compressed.mat"))
##'
##' unlink("test-uncompressed.mat")
##' unlink("test-compressed.mat")
##'
##' ## Example how to read and write a S4 class with rmatio
##' ## Create 'DemoS4Mat' class
##' setClass("DemoS4Mat",
##'          representation(a = "dgCMatrix",
##'                         b = "integer",
##'                         c = "matrix",
##'                         d = "numeric"))
##'
##' ## Create a function to coerce a 'DemoS4Mat' object to a list.
##' setAs(from="DemoS4Mat",
##'       to="list",
##'       def=function(from)
##'       {
##'         return(list(a=from@@a,
##'                     b=from@@b,
##'                     c=from@@c,
##'                     d=from@@d))
##'       }
##' )
##'
##' ## Create a function to coerce a list to a 'DemoS4Mat' object.
##' setAs(from="list",
##'       to="DemoS4Mat",
##'       def=function(from)
##'       {
##'         return(new("DemoS4Mat",
##'                     a=from[["a"]],
##'                     b=as.integer(from[["b"]]),
##'                     c=from[["c"]],
##'                     d=from[["d"]]))
##'       }
##' )
##'
##' ## Define a method to write a 'DemoS4Mat' object to a MAT file.
##' setMethod("write.mat",
##'           signature(object = "DemoS4Mat"),
##'           function(object,
##'                    filename,
##'                    compression,
##'                    version)
##'           {
##'             ## Coerce the 'DemoS4Mat' object to a list and
##'             ## call 'rmatio' 'write.mat' with the list.
##'             return(write.mat(as(object, "list"),
##'                              filename,
##'                              compression,
##'                              version))
##'           }
##' )
##'
##' ## Create a new 'DemoS4Mat' object
##' demoS4mat <- new("DemoS4Mat",
##'                  a = Matrix(c(0, 0, 0, 0, 0, 0, 1, 0, 0,
##'                               0, 0, 0, 0, 0, 0, 0, 1, 0,
##'                               0, 0, 0, 0, 0, 0, 0, 0, 1),
##'                               nrow=3,
##'                               ncol=9,
##'                               byrow=TRUE,
##'                               sparse=TRUE),
##'                  b = 1:5,
##'                  c = matrix(as.numeric(1:9), nrow=3),
##'                  d = c(6.0, 7.0, 8.0))
##'
##' ## Write to MAT file
##' write.mat(demoS4mat, filename)
##'
##' ## Read the MAT file
##' demoS4mat.2 <- as(read.mat(filename), "DemoS4Mat")
##'
##' ## Check result
##' stopifnot(identical(demoS4mat, demoS4mat.2))
##'
##' unlink(filename)
##' }
setGeneric("write.mat",
           signature = "object",
           function(object,
                    filename = NULL,
                    compression = TRUE,
                    version = c('MAT5')) standardGeneric("write.mat"))

##' @rdname write.mat-methods
##' @include have_lib.r
##' @export
setMethod("write.mat",
          signature(object = "list"),
          function(object,
                   filename,
                   compression,
                   version)
          {
            ## Check filename
            if(any(!is.character(filename),
                   !identical(length(filename), 1L),
                   nchar(filename) < 1)) {
              stop("'filename' must be a character vector of length one")
            }

            ## Check compression
            if(any(!is.logical(compression),
                   !identical(length(compression), 1L))) {
              stop("'compression' must be a logical vector of length one")
            }

            if(identical(compression, TRUE)) {
                if(!have.zlib()) {
                    stop(paste("Sorry, library 'zlib' is not available.",
                               "Use 'compression=FALSE' or install with 'zlib'"))
                }
                compression = 1L
            } else {
                compression = 0L
            }

            ## Check version
            version <- match.arg(version)
            if(identical(version, 'MAT5')) {
              version <- 0x0100L
              header <- sprintf("MATLAB 5.0 MAT-file, Platform: %s, Created By: rmatio v%s on %s",
                                R.version$platform[[1]],
                                packageVersion('rmatio'),
                                date())
            } else {
              stop('Unsupported version')
            }

            ## Check names in object
            if(any(is.null(names(object)),
                   !all(nchar(names(object))),
                   any(duplicated(names(object))))) {
              stop("All values in the list must have a unique name")
            }

            .Call("write_mat", object, filename, compression, version, header)

            invisible(NULL)
          }
)
