##
##
## Copyright (c) 2009, Brandon Whitcher and Volker Schmid
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## $Id: dwi.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## setGeneric("ADC.fast")
#############################################################################
#' Estimate the Apparent Diffusion Coefficient (ADC)
#' 
#' Estimation of apparent diffusion coefficient (ADC) values, using a single
#' exponential function, is achieved through nonlinear optimization.
#' 
#' The \code{adc.lm} function estimates parameters for a vector of observed MR
#' signal intensities using the following relationship \deqn{S(b) = S_0
#' \exp(-bD),} where \eqn{S_0}{S0} is the baseline signal intensity and \eqn{D}
#' is the apparent diffusion coefficient (ADC).  It requires the routine
#' \code{nls.lm} that applies the Levenberg-Marquardt algorithm.  Note, low
#' b-values (\eqn{<50} or \eqn{<100} depending on who you read) should be
#' avoided in the parameter estimation because they do not represent
#' information about the diffusion of water in tissue.
#' 
#' The \code{ADC.fast} function rearranges the assumed multidimensional (2D or
#' 3D) structure of the DWI data into a single matrix to take advantage of
#' internal R functions instead of loops, and called \code{adc.lm}.
#' 
#' @aliases adc.lm ADC.fast,array-method ADC.fast
#' @param signal Signal intensity vector as a function of b-values.
#' @param b,bvalues Diffusion weightings (b-values).
#' @param guess Initial values of \eqn{S_0}{S0} and \eqn{D}.
#' @param control An optional list of control settings for \code{nls.lm}.  See
#' \code{nls.lm.control} for the names of the settable control values and their
#' effect.
#' @param dwi Multidimensional array of diffusion-weighted images.
#' @param ... Additional variables defined by the method.  
#' @param dwi.mask Logical array that defines the voxels to be analyzed.
#' @param multicore is a logical variable (default = \code{FALSE}) that allows
#' parallel processing via \pkg{parallel}.
#' @param verbose Additional information will be printed when
#' \code{verbose=TRUE}.
#' @return A list structure is produced with estimates of \eqn{S_0}, \eqn{D}
#' and information about the convergence of the nonlinear optimization routine.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link[minpack.lm]{nls.lm}}
#' @references Buxton, R.B. (2002) \emph{Introduction to Functional Magnetic
#' Resonance Imaging: Principles & Techniques}, Cambridge University Press:
#' Cambridge, UK.
#' 
#' Callahan, P.T. (2006) \emph{Principles of Nuclear Magnetic Resonance
#' Microscopy}, Clarendon Press: Oxford, UK.
#' 
#' Koh, D.-M. and Collins, D.J. (2007) Diffusion-Weighted MRI in the Body:
#' Applications and Challenges in Oncology, \emph{American Journal of
#' Roentgenology}, \bold{188}, 1622-1635.
#' @keywords models
#' @examples
#' 
#' S0 <- 10
#' b <- c(0, 50, 400, 800)  # units?
#' D <- 0.7e-3              # mm^2 / sec (normal white matter)
#' 
#' ## Signal intensities based on the (simplified) Bloch-Torry equation
#' dwi <- function(S0, b, D) {
#'   S0 * exp(-b*D)
#' }
#' 
#' set.seed(1234)
#' signal <- array(dwi(S0, b, D) + rnorm(length(b), sd=0.15),
#'                 c(rep(1,3), length(b)))
#' ADC <- ADC.fast(signal, b, array(TRUE, rep(1,3)))
#' unlist(ADC) # text output
#' 
#' par(mfrow=c(1,1)) # graphical output
#' plot(b, signal, xlab="b-value", ylab="Signal intensity")
#' lines(seq(0,800,10), dwi(S0, seq(0,800,10), D), lwd=2, col=1)
#' lines(seq(0,800,10), dwi(c(ADC$S0), seq(0,800,10), c(ADC$D)), lwd=2, col=2)
#' legend("topright", c("True","Estimated"), lwd=2, col=1:2)
#' 
#' @export
#' @docType methods
#' @rdname ADC-methods
setGeneric("ADC.fast", function(dwi, ...) standardGeneric("ADC.fast"))
#' @export
#' @rdname ADC-methods
#' @aliases ADC.fast,array-method
setMethod("ADC.fast", signature(dwi="array"),
          function(dwi, bvalues, dwi.mask,
                   control=minpack.lm::nls.lm.control(maxiter=150),
                   multicore=FALSE, verbose=FALSE)
          .dcemriWrapper("ADC.fast", dwi, bvalues, dwi.mask, control,
                         multicore, verbose))

.ADC.fast <- function(dwi, bvalues, dwi.mask,
                      control=minpack.lm::nls.lm.control(maxiter=150),
                      multicore=FALSE, verbose=FALSE) {
  if (length(dim(dwi)) != 4) { # Check dwi is a 4D array
    stop("Diffusion-weighted data must be a 4D array.")
  }
  if (!is.logical(dwi.mask)) { # Check dyn.mask is logical
    stop("Mask must be logical.")
  }
  nvalues <- length(bvalues)
  nvoxels <- sum(dwi.mask)
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  dwi.mat <- matrix(dwi[dwi.mask], nvoxels)
  dwi.list <- vector("list", nvoxels)
  for (k in 1:nvoxels) {
    dwi.list[[k]] <- dwi.mat[k,]
  }
  if (verbose) {
    cat("  Calculating S0 and D...", fill=TRUE)
  }
  if (multicore) {
    fit.list <- parallel::mclapply(dwi.list, function(x) {
      adc.lm(x, bvalues, guess=c(0.75*x[1], 0.001), control)
    })
  } else {
    fit.list <- lapply(dwi.list, function(x) {
      adc.lm(x, bvalues, guess=c(0.75*x[1], 0.001), control)
    })
  }
  rm(dwi.list) ; gc()
  S0 <- D <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  for (k in 1:nvoxels) {
    if (fit.list[[k]]$info > 0 && fit.list[[k]]$info < 5) {
      S0$par[k] <- fit.list[[k]]$S0
      D$par[k] <- fit.list[[k]]$D
      S0$error[k] <- sqrt(fit.list[[k]]$hessian[1,1])
      D$error[k] <- sqrt(fit.list[[k]]$hessian[2,2])
    } else {
      S0$par[k] <- D$par[k] <- S0$error[k] <- D$error[k] <- NA
    }
  }
  rm(fit.list) ; gc()
  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  S0.array <- D.array <- S0error <- Derror <- array(NA, dim(dwi)[1:3])
  S0.array[dwi.mask] <- S0$par
  D.array[dwi.mask] <- D$par
  S0error[dwi.mask] <- S0$error
  Derror[dwi.mask] <- D$error

  list(S0 = S0.array, D = D.array, S0.error = S0error, D.error = Derror)
}

#############################################################################
## adc.lm() = estimate ADC using Levenburg-Marquardt
#############################################################################
#' @rdname ADC-methods
#' @export
adc.lm <- function(signal, b, guess, control=minpack.lm::nls.lm.control()) {
  func <- function(x, y) {
    S0 <- x[1]
    D <- x[2]
    signal <- y[[1]]
    b <- y[[2]]
    signal - S0 * exp(-b*D)
  }
  out <- minpack.lm::nls.lm(par=guess, fn=func, control=control,
                            y=list(signal, b))
  list(S0=out$par[1], D=out$par[2], hessian=out$hessian, info=out$info,
       message=out$message)
}
