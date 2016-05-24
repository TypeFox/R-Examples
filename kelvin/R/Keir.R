#' Complementary solution to the Kelvin differential equation (K)
#'
#' This function calculates the complex solution to the Kelvin differential
#' equation using modified Bessel functions of the \emph{second kind}, specifically
#' those produced by \code{\link[Bessel]{BesselK}}.
#'
#' @details
#' \code{\link{Ker}} and \code{\link{Kei}} are wrapper functions
#' which return the real and imaginary components of \code{\link{Keir}},, respectively.
#' 
#' @param x numeric; values to evaluate the complex solution at
#' @param nu. numeric; value of \eqn{\nu} in \eqn{\mathcal{K}_\nu}{Kei,Ker} solutions
#' @param nSeq. positive integer; equivalent to \code{nSeq} in \code{\link[Bessel]{BesselK}}
#' @param add.tol logical; Should a fudge factor be added to prevent an error for zero-values?
#' @param return.list logical; Should the result be a list instead of matrix?
#' @param show.scaling logical; Should the normalization values be given as a message?
#' @param ... additional arguments passed to \code{\link[Bessel]{BesselK}} or \code{\link{Keir}}
#' 
#' @export
#' @name Keir
#' 
#' @return If \code{return.list==FALSE} (the default),
#' a complex matrix with as many columns as using \code{nSeq.} creates.
#' Otherwise the result is a list with matrices for
#' Real and Imaginary components.
#'
#' @author Andrew Barbour
#' 
#' @references \url{http://mathworld.wolfram.com/KelvinFunctions.html}
#' @references Imaginary: \url{http://mathworld.wolfram.com/Kei.html}
#' @references Real: \url{http://mathworld.wolfram.com/Ker.html}
#' 
#' @seealso \code{\link{kelvin-package}}, \code{\link{Beir}}, \code{\link[Bessel]{BesselK}}
#' 
#' @examples
#' 
#' Keir(1:10)    # defaults to nu.=0, nSeq=1
#' Keir(1:10, nu.=2)
#' Keir(1:10, nSeq=2)
#' Keir(1:10, nSeq=2, return.list=TRUE)
#' 
Keir <- function(x, ...) UseMethod("Keir")

#' @rdname Keir
#' @export
Keir.default <- function(x, nu.=0, nSeq.=1, add.tol=TRUE, return.list=FALSE, show.scaling=FALSE, ...){
  if (add.tol){
    ret.ind <- FALSE
    #heuristic fix for zero values
    tol <- 1e-12
    zero.inds <- x == 0
    if (any(zero.inds)){
      ret.ind <- TRUE
      warning(sprintf('values of zero were replaced with  %s',tol))
      x[zero.inds] <- tol
    }
  } else {
    if (0 %in% x) stop("zeros in 'x'")
  }
  #
  BessX <- x * exp(pi * complex(real=0, imaginary = 1) / 4)
  # Bug fix: must multiply by the specific scaling for nu., so if
  # nSeq is given the scaling will be wrong.  Fix is to create a
  # vector of scalings.  This page was useful:
  #http://keisan.casio.com/has10/SpecExec.cgi
  Nu. <- nu.:(nu.+nSeq.-1)
  Bsc <- zapsmall(exp(-1 * pi * Nu. * complex(real=0, imaginary = 1) / 2))
  if (show.scaling) {message(sprintf("\t>>>>\tnu=%i\tscaling:\t%s\n", Nu., Bsc))}
  #
  Bsl <- Bessel::BesselK(BessX, nu=nu., nSeq=nSeq., ...)
  nr. <- length(as.vector(BessX))
  stopifnot(!is.null(nr.))
  #
  Bsl <- Bsl*matrix(rep(Bsc, nr.), nrow=nr., byrow=TRUE)
  #
  if (return.list){
    Bsl <- list(kei=Im(Bsl), ker=Re(Bsl))
    if (ret.ind) Bsl[['zero.indices']] <- zero.inds
  }
  return(Bsl)
}

#' @rdname Keir
#' @export
#' @examples
#' 
#' # Imaginary component only
#' Kei(1:10)
Kei <- function(...) Keir(..., return.list=TRUE)[['kei']]

#' @rdname Keir
#' @export
#' @examples
#' 
#' # Real component only
#' Ker(1:10)
Ker <- function(...) Keir(..., return.list=TRUE)[['ker']]
