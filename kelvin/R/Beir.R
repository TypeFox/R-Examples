#' Fundamental solution to the Kelvin differential equation (J)
#'
#' This function calculates the complex solution to the Kelvin differential
#' equation using modified Bessel functions of the \emph{first kind}, specifically
#' those produced by \code{\link[Bessel]{BesselJ}}.
#'
#' @details
#' \code{\link{Ber}} and \code{\link{Bei}} are wrapper functions
#' which return the real and imaginary components of \code{\link{Beir}}, respectively.
#' 
#' @inheritParams Keir
#' @param nu. numeric; value of \eqn{\nu} in \eqn{\mathcal{B}_\nu}{Bei,Ber} solutions
#' @param nSeq. positive integer; equivalent to \code{nSeq} in \code{\link[Bessel]{BesselJ}}
#' @param ... additional arguments passed to \code{\link[Bessel]{BesselK}} or \code{\link{Beir}}
#' 
#' @export
#' @name Beir
#'
#' @return If \code{return.list==FALSE} (the default),
#' a complex matrix with as many columns as using \code{nSeq.} creates.
#' Otherwise the result is a list with matrices for
#' Real and Imaginary components.
#' 
#' @author Andrew Barbour
#' 
#' @references \url{http://mathworld.wolfram.com/KelvinFunctions.html}
#' @references Imaginary: \url{http://mathworld.wolfram.com/Bei.html}
#' @references Real: \url{http://mathworld.wolfram.com/Ber.html}
#' 
#' @seealso \code{\link{kelvin-package}}, \code{\link{Keir}}, \code{\link[Bessel]{BesselJ}}
#' 
#' @examples
#' 
#' Beir(1:10)    # defaults to nu.=0
#' Beir(1:10, nu.=2)
#' Beir(1:10, nSeq.=2)
#' Beir(1:10, nSeq.=2, return.list=TRUE)
#' 
Beir <- function(x, ...) UseMethod("Beir")

#' @rdname Beir
#' @export
Beir.default <- function(x, nu.=0, nSeq.=1, return.list=FALSE, ...){
  bess <- Bessel::BesselJ(x * exp(3 * pi * complex(real=0, imaginary = 1) / 4), nu=nu., nSeq=nSeq., ...)
  bess <- as.matrix(bess)
  if (return.list){
    bess <- list(bei=Im(bess), ber=Re(bess))
  }
  return(bess)
}

#' @rdname Beir
#' @export
#' @examples
#' 
#' # Imaginary component only
#' Bei(1:10)
Bei <- function(...) Beir(..., return.list=TRUE)[['bei']]

#' @rdname Beir
#' @export
#' @examples
#' 
#' # Real component only
#' Ber(1:10)
Ber <- function(...) Beir(..., return.list=TRUE)[['ber']]
