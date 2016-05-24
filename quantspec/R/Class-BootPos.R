################################################################################
#' Class for Generation of Bootstrapped Replications of a Time Series.
#'
#' \code{BootPos} is an S4 class that provides a common interface
#' to different algorithms that can be used for implementation of a block
#' bootstrap procedure in the time domain.
#'
#' After initialization the bootstrapping can be performed by applying
#' \code{getPositions} to the object.
#'
#' Different block bootstraps are implemented by creating a subclass together
#' with a \code{getPositions} method that contains the implementation of the
#' block resampling procedure.
#'
#' Currently the following implementations are available:
#'
#' \itemize{
#' 		\item \code{\link{MovingBlocks}} and \code{\link{getPositions-MovingBlocks}}.
#' }
#'
#' @name   BootPos-class
#' @aliases BootPos
#' @exportClass BootPos
#'
#' @keywords S4-classes
#'
#' @slot l the (expected) block length for the block bootstrap methods
#' @slot N number of available observations to bootstrap from
#'
#' @references
#' Lahiri, S. N. (1999). Theoretical Comparisons of Block Bootstrap Methods.
#' \emph{The Annals of Statistics}, \bold{27}(1), 386--404.
################################################################################

setClass(
    Class = "BootPos",
    representation=representation(
        l = "numeric",
        N = "numeric"
    )
)
