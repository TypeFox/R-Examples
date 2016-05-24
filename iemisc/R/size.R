#' Size of R objects (GNU Octave/MATLAB compatible)
#'
#' Provides the dimensions of R objects in a manner compatible with
#' GNU Octave/MATLAB. This function is the same as \code{\link[pracma]{size}}, except this
#' \code{size} can find the size of character vectors too. Some documentation from
#' \code{\link[pracma]{size}}.
#'
#'
#' @param x An R object (array, vector, or matrix)
#' @param k integer specifying a particular dimension
#'
#' @return "Return the number of rows and columns of the object x as a numeric
#'   vector. If given a second argument, \code{size} will return the size of the
#'   corresponding dimension." Source: Eaton.
#'
#'
#' @source
#' pracma size function definition - R package pracma created and maintained by Hans Werner Borchers. See \code{\link[pracma]{interp1}}.
#'
#'
#' @references
#' John W. Eaton, David Bateman, and Søren Hauberg (2009). \emph{GNU Octave version 3.0.1 manual: a high-level interactive language for numerical computations}. CreateSpace Independent Publishing Platform. ISBN 1441413006, URL \url{http://www.gnu.org/software/octave/doc/interpreter/}. Page 42.
#'
#'
#'
#' @author Hans Werner Borchers (pracma size), Irucka Embry
#'
#'
#'
#' @encoding UTF-8
#'
#'
#'
#'
#' @seealso \code{\link[base]{dim}}, \code{\link[pracma]{size}}
#'
#'
#' @examples
#' library(iemisc)
#' # Examples from GNU Octave size
#' object1 <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2, byrow = TRUE)
#' size(object1)
#'
#' source("http://gsubfn.googlecode.com/svn/trunk/R/list.R")
#' list[nr, nc] <- size(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2,
#'                 byrow = TRUE))
#'
#' size(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2, byrow = TRUE), 2)
#'
#' # Examples from pracma size
#' size(1:8)
#'
#' size(matrix(1:8, 2, 4))
#'
#' size(matrix(1:8, 2, 4), 2)
#'
#' size(matrix(1:8, 2, 4), 3)
#'
#' ss <- "object"
#' size(ss)
#'
#'
#' \dontrun{
#' # check against GNU Octave
#' library(RcppOctave) # requires Octave (>= 3.2.4) and its development files
#' o_source(text = "
#' \% Examples from GNU Octave size
#' object1 = [1, 2; 3, 4; 5, 6];
#' size(object1)
#'
#' [nr, nc] = size([1, 2; 3, 4; 5, 6])
#'
#' size([1, 2; 3, 4; 5, 6], 2)
#'
#' \% Examples from pracma size
#' size(1:8)
#'
#' object2 = [1 3 5 7; 2 4 6 8];
#'
#' size(object2)
#'
#' size(object2, 2)
#'
#' size(object2, 3)
#'
#' ss = 'object';
#' size(ss)
#' ")
#' }
#'
#'
#' @export
# Sources 1 begins
size <- function (x, k)
{
    if (length(x) == 0)

        sz <- 0

    else if (is.vector(x) & is.vector(x, mode = "character") == FALSE)

    sz <- c(1, length(x))

	else if (is.vector(x) & is.vector(x, mode = "character"))

	sz <- c(1, nchar(x))

    else if (is.array(x))

    sz <- dim(x)

    else sz <- NULL

    if (!missing(k)) {

    if (k > length(sz))

    sz <- 1

    else if (k >= 1)

    sz <- sz[k]

    else stop("Requested dimension 'k' is out of range.")

    }

    return(sz)
# Sources 1 ends
    }



#' Length of R objects (GNU Octave/MATLAB compatible)
#'
#' Obtain the length of R objects [arrays, matrices, and vectors (including
#' lists)] in a manner compatible with GNU Octave/MATLAB. Some documentation
#' from \code{\link[base]{length}}.
#'
#'
#' @param x An R object (array, matrix, vector)
#'
#' @return Return the length of the object x as an integer. "The length is 0
#'   for empty objects, 1 for scalars (in R, a \code{vector} of \code{length} 1), and
#'   the number of elements (in R, the \code{length}) for vectors. For matrix objects,
#'   the length is the number of rows or columns, whichever is greater (this
#'   odd definition is used for compatibility with MATLAB)." Source: Eaton.
#'
#' @references
#' \enumerate{
#'    \item Samit Basu (2002-2006). FreeMat v4.0, \url{http://freemat.sourceforge.net/}.
#'    \item John W. Eaton, David Bateman, and Søren Hauberg (2009). \emph{GNU Octave version 3.0.1 manual: a high-level interactive language for numerical computations}. CreateSpace Independent Publishing Platform. ISBN 1441413006, URL \url{http://www.gnu.org/software/octave/doc/interpreter/}. Page 41.
#' }
#'
#'
#' @author Irucka Embry, Samit Basu (FreeMat)
#'
#'
#' @encoding UTF-8
#'
#'
#'
#'
#' @seealso \code{\link[base]{length}}, \code{\link[pracma]{size}}, \code{\link{size}}
#'
#'
#' @examples
#' library(iemisc)
#' library(import)
#' import::from(pracma, ones)
#' # Example from pracma isempty
#' object1 <- matrix(0, 1, 0)
#' lengths(object1)
#'
#' object2 <- 2
#' lengths(object2)
#'
#' object3 <- 1:10
#' lengths(object3)
#'
#' object4 <- ones(3, 4)
#' lengths(object4)
#'
#' object5 <- "ss"
#' lengths(object5)
#'
#' object6 <- list(letters, b <- 2)
#' lengths(object6)
#'
#'
#' \dontrun{
#' # check against GNU Octave
#' library(RcppOctave) # requires Octave (>= 3.2.4) and its development files
#' o_source(text = "
#' object1 = [];
#' length(object1)
#'
#' object2 = 2;
#' length(object2)
#'
#' object3 = 1:10;
#' length(object3)
#'
#' object4 = ones(3, 4);
#' length(object4)
#'
#' object5 = 'ss';
#' length(object5)
#' ")
#' }
#'
#'
#' @importFrom pracma isempty
#'
#' @export
lengths <- function (x) {
if (isempty(x))

  0

else if (length(is.vector(x)) == 1)

  1

else if (is.vector(x))

  length(x)

else (is.matrix(x))

  max(size(x))

  }



#' Number of elements (GNU Octave/MATLAB compatible)
#'
#' Obtain the number of elements of R objects [arrays, matrices, and vectors
#' (including lists)] in a manner compatible with GNU Octave/MATLAB. Some
#' documentation from \code{\link[base]{length}}.
#'
#'
#' @param x An R object (array, matrix, vector)
#' @param ... R objects (indices idx1, idx2, ...)
#'
#' @return "Return the number of elements in the R object x. Optionally, if
#'   indices idx1, idx2, ... are supplied, return the number of elements that
#'   would result from the indexing a(idx1, idx2, ...)." Source: Eaton page 41.
#'
#'
#' @source
#' \enumerate{
#'    \item r - Add a Column to a Dataframe From a List of Values - Stack Overflow answered by Matthew Plourde on Jun 21 2012. See \url{http://stackoverflow.com/questions/11130037/add-a-column-to-a-dataframe-from-a-list-of-values/11130178}.
#'    \item r - Why does is.vector() return TRUE for list? - Stack Overflow answered by Andrie on May 17 2011. See \url{http://stackoverflow.com/questions/6032772/why-does-is-vector-return-true-for-list/6032909}.
#' }
#'
#'
#' @references
#' \enumerate{
#'    \item Samit Basu (2002-2006). FreeMat v4.0, \url{http://freemat.sourceforge.net/}.
#'    \item John W. Eaton, David Bateman, and Søren Hauberg (2009). \emph{GNU Octave version 3.0.1 manual: a high-level interactive language for numerical computations}. CreateSpace Independent Publishing Platform. ISBN 1441413006, URL \url{http://www.gnu.org/software/octave/doc/interpreter/}. Page 41.
#' }
#'
#'
#' @author Irucka Embry, Samit Basu (FreeMat)
#'
#'
#' @encoding UTF-8
#'
#'
#'
#'
#' @seealso \code{\link[matlab]{numel}}, \code{\link[pracma]{numel}}, \code{\link{size}}, \code{\link{length}}
#'
#'
#' @examples
#' library(iemisc)
#' library(import)
#' import::from(pracma, ones)
#' xx <- list(1:26, 1:10)
#' numel(xx)
#'
#' # Examples from GNU Octave numel
#' a <- 1
#' b <- ones(2, 3)
#' numel(a, b)
#'
#' a <- 2
#' b <- ones(2, 3)
#' c <- ones(3, 4)
#' numel(a, b)
#' numel(a, b, c)
#'
#' f <- matrix(c(10, 12, 23, 21, 62, 93), nrow = 2, ncol = 3, byrow = TRUE)
#' g <- c(2, 4)
#' numel(f, g)
#'
#'
#' \dontrun{
#' # check against GNU Octave
#' library(RcppOctave) # requires Octave (>= 3.2.4) and its development files
#' o_source(text = "
#' xx = {1:26, 1:10}
#'
#' \% Examples from GNU Octave numel
#' a = 1;
#' b = ones(2, 3);
#' numel(a, b)
#'
#' a = 2;
#' b = ones(2, 3);
#' c = ones(3, 4);
#' numel(a, b)
#' numel(a, b, c)
#'
#' f = [10 12 23; 21 62 93];
#' g = [2 4];
#' numel(f, g)
#' ")
#' }
#'
#'
#' @export
numel <- function (x, ...) {

if (nargs() == 1) {

  lens <- prod(size(x))

} else {

  varargin <- list(...)

  len <- 1

  lens <- vector("list", numel(varargin))
# Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (k in 1:numel(varargin)) {

  lens[[k]] <- len * size(varargin[[k]])
}

  lens <- prod(unlist(lens))

}

return(lens)

}
