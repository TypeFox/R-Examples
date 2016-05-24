#' @include Class-BootPos.R
NULL

################################################################################
#' Class for Moving Blocks Bootstrap implementation.
#'
#' \code{MovingBlocks} is an S4 class that implements the moving blocks
#' bootstrap described in K{\"u}nsch (1989).
#'
#' \code{MovingBlocks} extends the S4 class
#' \code{\link{BootPos}} and the remarks made in its documentation
#' apply here as well.
#'
#' The Moving Blocks Bootstrap method of K{\"u}nsch (1989) resamples blocks
#' randomly, with replacement from the collection of overlapping blocks of
#' length \code{l} that start with observation 1, 2, \ldots, \code{N-l+1}.
#' A more precise description of the procedure can also be found in
#' Lahiri (1999), p. 389.
#'
#' @name   MovingBlocks-class
#' @aliases MovingBlocks
#' @exportClass MovingBlocks
#'
#' @encoding latin1
#'
#' @keywords S4-classes
#'
#' @seealso \code{\link{getPositions-MovingBlocks}}
#'
#' @references
#' K{\"u}nsch, H. R. (1989). The jackknife and the bootstrap for general stationary
#' observations. \emph{The Annals of Statistics}, \bold{17}, 1217--1261.
################################################################################

setClass(
    Class = "MovingBlocks",
    contains = "BootPos"
)

setMethod(
    f = "initialize",
    signature = "MovingBlocks",
    definition = function(.Object, l, N) {

      .Object@l <- l
      .Object@N <- N

      # Return object
      return(.Object)
    }
)

################################################################################
#' Get Positions for the Moving Blocks Bootstrap.
#'
#' @name getPositions-MovingBlocks
#' @aliases getPositions,MovingBlocks-method
#' 
#' @importFrom stats runif
#'
#' @param object a \code{MovingBlocks} object; used to specify the parameters
#'                \code{N}, \code{l} and the type of the bootstrap.
#' @param B Number of independent repetitions to bootstrap.
#'
#' @return a matrix of dimension \code{[N,B]} where each column gives the
#'         positions in which to reorder the observations to yield one
#'          bootstrap replication.
################################################################################

setMethod(f = "getPositions",
    signature = "MovingBlocks",
    definition = function(object, B=1) {

    N <- object@N
    l <- object@l
    nBlocks <- ceiling(N/l)

    positions <- c()

    for (b in 1:B) {
      blocks <- matrix(ncol=nBlocks, nrow=l)
      blocks[1,] <- floor(runif(n=nBlocks, min=1,max=N-l+1))
      if (l > 1) {
        for (i in 2:l) {
          blocks[i,] <- blocks[1,]+i-1
        }
      }
      positions <- c(positions,as.vector(blocks)[1:N])
    }

    return(matrix(positions,nrow=N))
  }
)


################################################################################
#' Create an instance of the \code{\link{MovingBlocks}} class.
#'
#' @name MovingBlocks-constructor
#' @aliases movingBlocks
#' @export
#'
#' @keywords Constructors
#'
#' @param l the block length for the block bootstrap methods
#' @param N number of available observations to bootstrap from
#'
#' @return Returns an instance of \code{MovingBlocks}.
################################################################################

movingBlocks <- function( l, N ) {

  if (!(is.wholenumber(l) && is.wholenumber(N) && 0 < l && l <= N)) {
    stop("'l' and 'N' need to be specified as integers with 0 < l <= N")
  }

  obj <- new(
      Class = "MovingBlocks",
      l = l,
      N = N
  )

  return(obj)
}
