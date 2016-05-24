##
## Copyright (c) 2009,2010 Brandon Whitcher and Volker Schmid
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
## $Id:$
##

#############################################################################
## setGeneric("T2.fast")
#############################################################################
#' Quantitative T2 Methods
#' 
#' The regional blood volume is found by integrating of the tissue concentration
#' curve and the artieral input funciton (AIF).  In order to avoid reperfusion
#' effects on the rCBV measurements, the tissue and arteiral concentration 
#' curves must first be reduced to their first-pass versions.  
#' 
#' @aliases T2.lm T2.fast T2.fast,array-method T2.fast,anlz-method
#' T2.fast,nifti-method
#' @param signal is the vector of signal intensities as a function of echo
#' times.
#' @param TE is the vector of echo times (in seconds).
#' @param guess is the vector of initial values for the parameters of interest:
#' \eqn{\rho}{rho} and \eqn{T2}{T2}.
#' @param control An optional list of control settings for \code{nls.lm}.  See
#' \code{nls.lm.control} for the names of the settable control values and their
#' effect.
#' @param cpmg is a multidimensional array of signal intensities.  The last
#' dimension is assumed to be a function of the echo times, while the previous
#' dimenions are assued to be spatial.
#' @param ... Additional variables defined by the method.  
#' @param cpmg.mask is a (logical) multidimensional array that identifies the
#' voxels to be analyzed.
#' @param multicore is a logical variable (default = \code{FALSE}) that allows
#' parallel processing via \pkg{multicore}.
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @return A list structure is produced with (all or some of the) parameter
#' estimates 
#' \item{rho}{Scaling factor between signal intensity and T2 (proton density).} 
#' \item{T2}{T2 relaxation time.}
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{R1.fast}}, \code{\link{R10.lm}}
#' @references 
#' Kennan, R.P. and J\"ager, H.R. (2004) $T_2$- and $T_2^*$-w DCE-MRI: Blood
#' Perfusion and Volume Estimation using Bolus Tracking, in \emph{Quantiative 
#' MRI of the Brain} (P. Tofts ed.), Wiley: Chichester, UK, pp. 365-412.
#' @keywords misc
#' @rdname CPMG-methods
#' @export
#' @docType methods
setGeneric("T2.fast", function(cpmg, ...) standardGeneric("T2.fast"))
#' @export
#' @rdname CPMG-methods
#' @aliases T2.fast,array-method
setMethod("T2.fast", signature(cpmg="array"),
          function(cpmg, cpmg.mask, TE,
                   control=minpack.lm::nls.lm.control(maxiter=150),
                   multicore=FALSE, verbose=FALSE)
          .dcemriWrapper("T2.fast", cpmg, cpmg.mask, TE, control, multicore,
                         verbose))

#############################################################################
## T2.fast()
#############################################################################

.T2.fast <- function(cpmg, cpmg.mask, TE,
                     control=minpack.lm::nls.lm.control(maxiter=150),
                     multicore=FALSE, verbose=FALSE) {

  if (length(dim(cpmg)) != 4) { # Check cpmg is a 4D array
    stop("CPMG data must be a 4D array.")
  }
  if (!is.logical(cpmg.mask)) { # Check cpmg.mask is logical
    stop("Mask must be logical.")
  }
  X <- nrow(cpmg)
  Y <- ncol(cpmg)
  Z <- oro.nifti::nsli(cpmg)
  nvoxels <- sum(cpmg.mask)
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  cpmg.mat <- matrix(cpmg[cpmg.mask], nrow=nvoxels)
  cpmg.list <- vector("list", nvoxels)
  for (k in 1:nvoxels) {
    cpmg.list[[k]] <- cpmg.mat[k,]
  }
  if (verbose) {
    cat("  Calculating T2 and rho...", fill=TRUE)
  }
  if (multicore) {
    T2.list <- parallel::mclapply(cpmg.list, function(x) {
      T2.lm(x, TE, guess=c(0.75 * x[1], 0.05), control)
    })
  } else {
    T2.list <- lapply(cpmg.list, function(x) {
      T2.lm(x, TE, guess=c(0.75 * x[1], 0.05), control)
    })
  }
  rm(cpmg.list) ; gc()
  T2 <- rho <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  pb <- txtProgressBar()
  for (k in 1:nvoxels) {
    if (T2.list[[k]]$info > 0 && T2.list[[k]]$info < 5) {
      T2$par[k] <- T2.list[[k]]$T2
      rho$par[k] <- T2.list[[k]]$rho
      T2$error[k] <- sqrt(T2.list[[k]]$hessian[1,1])
      rho$error[k] <- sqrt(T2.list[[k]]$hessian[2,2])
    } else {
      T2$par[k] <- rho$par[k] <- T2$error[k] <- rho$error[k] <- NA
    }
    setTxtProgressBar(pb, k)
  }
  close(pb)
  rm(T2.list) ; gc()
  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  T2.array <- rho.array <- T2.error <- rho.error <- array(NA, dim(cpmg)[1:3])
  T2.array[cpmg.mask] <- T2$par
  rho.array[cpmg.mask] <- rho$par
  T2.error[cpmg.mask] <- T2$error
  rho.error[cpmg.mask] <- rho$error

  list(rho = rho.array, T2 = T2.array, rho.error = rho.error,
       T2.error = T2.error)
}

#############################################################################
## T2.lm() = estimate exp(-TE/T2) using Levenburg-Marquardt
#############################################################################
#' @export
#' @rdname CPMG-methods
#' @aliases T2.lm
T2.lm <- function(signal, TE, guess, control=minpack.lm::nls.lm.control()) {
  func <- function(x, y) {
    rho <- x[1]
    T2 <- x[2]
    signal <- y[[1]]
    TE <- y[[2]]
    signal - rho * exp(-TE/T2)
  }
  out <- minpack.lm::nls.lm(par=guess, fn=func, control=control,
                            y=list(signal, TE))
  list(rho=out$par[1], T2=out$par[2], hessian=out$hessian, info=out$info,
       message=out$message)
}
