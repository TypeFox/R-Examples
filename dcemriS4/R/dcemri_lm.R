##
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
## $Id: dcemri_bayes.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## setGeneric("dcemri.lm")
#############################################################################
#' Pharmacokinetic Models for Dynamic Contrast-Enhanced MRI Data
#' 
#' Parameter estimation for single compartment models is performed using
#' literature-based or user-specified arterial input functions.  The
#' Levenburg-Marquardt algorithm does the heavy lifting.
#' 
#' Compartmental models are the solution to the modified general rate equation
#' (Kety 1951).  The specific parametric models considered here include the
#' basic Kety model
#' \deqn{C_t(t)=K^{trans}\left[C_p(t)\otimes\exp(-k_{ep}t)\right],} where
#' \eqn{\otimes} is the convoluation operator, and the so-called extended Kety
#' model
#' \deqn{C_t(t)=v_pC_p(t)+K^{trans}\left[C_p(t)\otimes\exp(-k_{ep}t)\right].}
#' The arterial input function must be either literature-based (with fixed
#' parameters) or the exponential AIF from Orton \emph{et al.} (2008) with
#' user-defined parameters.
#' 
#' @aliases dcemri.lm dcemri.lm,array-method
#' @param conc is a multidimensional (1D-4D) array of contrast agent
#' concentrations.  The last dimension is assumed to be temporal, while the
#' previous dimensions are assumed to be spatial.
#' @param time is a vector of acquisition times (in minutes) relative to
#' injection of the contrast agent.  Negative values should be used prior to
#' the injection.
#' @param img.mask is a (logical) multidimensional array that identifies the
#' voxels to be analyzed. Has to have same dimension as \code{conc} minus
#' temporal dimension.
#' @param model is a character string that identifies the type of compartmental
#' model to be used.  Acceptable models include: 
#' \describe{
#' \item{"weinmann"}{Tofts & Kermode AIF convolved with single
#' compartment model}
#' \item{"extended"}{Weinmann model extended with additional vascular 
#' compartment (default)}
#' \item{"orton.exp"}{Extended model using Orton's exponential AIF} 
#' \item{"orton.cos"}{Extended model using Orton's raised cosine AIF} 
#' \item{"kety.orton.exp"}{Kety model using Orton's exponential AIF} 
#' \item{"kety.orton.cos"}{Kety model using Orton's raised cosine AIF}
#' }
#' @param aif is a character string that identifies the parameters of the type
#' of arterial input function (AIF) used with the above model.  Acceptable
#' values are: 
#' \itemize{ 
#' \item\code{tofts.kermode}(default) for the \code{weinmann} and 
#' \code{extended} models 
#' \item\code{fritz.hansen} for the \code{weinmann} and \code{extended} models 
#' \item\dQuote{empirical} for the \code{weinmann} and \code{extended} models 
#' \item\code{orton.exp}(default) for the \code{orton.exp} and 
#' \code{kety.orton.exp} model
#' \item\code{orton.cos}(default) for the \code{orton.cos} and
#' \code{kety.orton.cos} model.  
#' \item\code{user} for the \code{orton.exp} and \code{orton.cos} model.
#' } 
#' All AIF models set the parametric form and parameter values -- except
#' \code{user}, where a set of user-defined parameter values are allowed, and
#' \code{empirical}, where a vector of values that fully characterize the
#' empirical AIF.
#' @param control is a list of parameters used by \code{nls.lm.control} that
#' are set by default, but may be customized by the user.
#' @param user is a list with the following parameters required: D, AB, muB,
#' AG, muG.
#' @param guess is a vector of starting values for kinetic parameter
#' estimation.  The vector must have length = 3 (with names \code{th0},
#' \code{th1} and \code{th3}) when the extended Kety model is used with the
#' vascular parameter and length = 2 (with names \code{th1} and \code{th3})
#' otherwise.
#' @param multicore is a logical variable (default = \code{FALSE}) that allows
#' parallel processing via \pkg{parallel}.
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @param ... Additional parameters to the function.
#' @return Parameter estimates and their standard errors are provided for the
#' masked region of the multidimensional array.  All multi-dimensional arrays
#' are provided in \code{nifti} format.
#' They include: 
#' \item{ktrans}{Transfer rate from plasma to the extracellular,
#' extravascular space (EES).} 
#' \item{kep}{Rate parameter for transport from the EES to plasma.} 
#' \item{ve}{Fractional occupancy by EES (the ratio between
#' \eqn{K^{trans}}{Ktrans} and \eqn{k_{ep}}{kep}).} 
#' \item{vp}{Fractional occupancy in the plasma space.} 
#' \item{ktranserror}{Standard error for \eqn{K^{trans}}{Ktrans}.} 
#' \item{keperror}{Standard error for \eqn{k_{ep}}{kep}.} 
#' \item{vperror}{Standard error for \eqn{v_p}{vp}.} 
#' The residual sum-of-squares is also provided, along with the original
#' acquisition times (for plotting purposes).
#' @note WARNING: when using the \code{empirical} AIF, a linear interpolation
#' is used to upsample the AIF to a one-second sampling rate.  This allows one
#' to utilize a computationally efficient numeric method to perform the
#' convolution.  If the empirical AIF is sampled faster than one Hertz, then
#' the upsampling operation will become a downsampling operation.  This should
#' not have any serious effect on the parameter estimates, but caution should
#' be exercised if very fast sampling rates are used to obtain an empirical
#' AIF.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com},\cr Volker
#' Schmid \email{volkerschmid@@users.sourceforge.net}
#' @seealso \code{\link{dcemri.bayes}}, \code{\link{dcemri.map}},
#' \code{\link{dcemri.spline}}, \code{\link[minpack.lm]{nls.lm}}
#' @references 
#' Ahearn, T.S., Staff, R.T., Redpath, T.W. and Semple, S.I.K.
#' (2005) The use of the Levenburg-Marquardt curve-fitting algorithm in
#' pharmacokinetic modelling of DCE-MRI data, \emph{Physics in Medicine and
#' Biology}, \bold{50}, N85-N92.
#' 
#' Fritz-Hansen, T., Rostrup, E., Larsson, H.B.W, Sondergaard, L., Ring, P. and
#' Henriksen, O. (1993) Measurement of the arterial concentration of Gd-DTPA
#' using MRI: A step toward quantitative perfusion imaging, \emph{Magnetic
#' Resonance in Medicine}, \bold{36}, 225-231.
#' 
#' Orton, M.R., Collins, D.J., Walker-Samuel, S., d'Arcy, J.A., Hawkes, D.J.,
#' Atkinson, D. and Leach, M.O. (2007) Bayesian estimation of pharmacokinetic
#' parameters for DCE-MRI with a robust treatment of enhancement onset time,
#' \emph{Physics in Medicine and Biology} \bold{52}, 2393-2408.
#' 
#' Orton, M.R., d'Arcy, J.A., Walker-Samuel, S., Hawkes, D.J., Atkinson, D.,
#' Collins, D.J. and Leach, M.O. (2008) Computationally efficient vascular
#' input function models for quantitative kinetic modelling using DCE-MRI,
#' \emph{Physics in Medicine and Biology} \bold{53}, 1225-1239.
#' 
#' Tofts, P.S., Brix, G, Buckley, D.L., Evelhoch, J.L., Henderson, E., Knopp,
#' M.V., Larsson, H.B.W., Lee, T.-Y., Mayr, N.A., Parker, G.J.M., Port, R.E.,
#' Taylor, J. and Weiskoff, R. (1999) Estimating kinetic parameters from
#' dynamic contrast-enhanced \eqn{T_1}-weighted MRI of a diffusable tracer:
#' Standardized quantities and symbols, \emph{Journal of Magnetic Resonance},
#' \bold{10}, 223-232.
#' 
#' Tofts, P.S. and Kermode, A.G. (1984) Measurement of the blood-brain barrier
#' permeability and leakage space using dynamic MR imaging. 1. Fundamental
#' concepts, \emph{Magnetic Resonance in Medicine}, \bold{17}, 357-367.
#' 
#' Weinmann, H.J., Laniado, M. and Mutzel, W. (1984) Pharmacokinetics of
#' Gd-DTPA/dimeglumine after intraveneous injection into healthy volunteers,
#' \emph{Physiological Chemistry and Physics and Medical NMR}, \bold{16},
#' 167-172.
#' @keywords models
#' @examples
#' 
#' data("buckley")
#' 
#' ## Empirical arterial input function
#' img <- array(t(breast$data), c(13,1,1,301))
#' time <- buckley$time.min
#' mask <- array(TRUE, dim(img)[1:3])
#' 
#' ## Estimate kinetic parameters directly from Buckley's empirical AIF
#' fit1 <- dcemri.lm(img, time, mask, model="weinmann", aif="empirical",
#'                   user=buckley$input)
#' fit2 <- dcemri.lm(img, time, mask, model="extended", aif="empirical",
#'                   user=buckley$input)
#' 
#' ## Set up breast data for dcemri
#' xi <- seq(5, 300, by=3)
#' img <- array(t(breast$data)[,xi], c(13,1,1,length(xi)))
#' time <- buckley$time.min[xi]
#' input <- buckley$input[xi]
#' 
#' ## Generate AIF params using the orton.exp function from Buckley's AIF
#' (aifparams <- orton.exp.lm(time, input))
#' fit3 <- dcemri.lm(img, time, mask, model="orton.exp", aif="user",
#'                   user=aifparams)
#' 
#' ## Scatterplot comparing true and estimated Ktrans values
#' plot(breast$ktrans, fit1$ktrans, xlim=c(0,0.75), ylim=c(0,0.75),
#'      xlab=expression(paste("True ", K^{trans})),
#'      ylab=expression(paste("Estimated ", K^{trans})))
#' points(breast$ktrans, fit2$ktrans, pch=2)
#' points(breast$ktrans, fit3$ktrans, pch=3)
#' abline(0, 1, lwd=1.5, col=2)
#' legend("bottomright", c("weinmann/empirical", "extended/empirical",
#'                         "orton.exp/user"), pch=1:3)
#' cbind(breast$ktrans, fit1$ktrans[,,1], fit2$ktrans[,,1], fit3$ktrans[,,1])
#' 
#' ## Scatterplot comparing true and estimated Ktrans values
#' plot(breast$vp, fit1$vp, type="n", xlim=c(0,0.15), ylim=c(0,0.15),
#'      xlab=expression(paste("True ", v[p])),
#'      ylab=expression(paste("Estimated ", v[p])))
#' points(breast$vp, fit2$vp, pch=2)
#' points(breast$vp, fit3$vp, pch=3)
#' abline(0, 1, lwd=1.5, col=2)
#' legend("bottomright", c("extended/empirical","orton.exp/user"), pch=2:3)
#' cbind(breast$vp, fit2$vp[,,1], fit3$vp[,,1])
#' 
#' @export
#' @docType methods
#' @import methods
#' @rdname dcemri.lm-methods
setGeneric("dcemri.lm", function(conc,  ...) standardGeneric("dcemri.lm"))
#' @export
#' @rdname dcemri.lm-methods
#' @aliases dcemri.lm,array-method
setMethod("dcemri.lm", signature(conc="array"),
	  function(conc,time,img.mask, model="extended", aif=NULL,
                   control=minpack.lm::nls.lm.control(), user=NULL,
                   guess=NULL, multicore=FALSE, verbose=FALSE, ...)
          .dcemriWrapper("dcemri.lm", conc, time, img.mask, model, aif,
                         control, user, guess, multicore, verbose, ...))

.dcemri.lm <- function(conc, time, img.mask, model="extended", aif=NULL,
                       control=minpack.lm::nls.lm.control(), user=NULL,
                       guess=NULL, multicore=FALSE, verbose=FALSE, ...) {
  switch(model,
         weinmann = ,
         extended = {
           aif <- ifelse(is.null(aif), "tofts.kermode", aif)
           aif.names <- c("tofts.kermode","fritz.hansen","empirical","user")
           if (! aif %in% aif.names) {
             stop(sprintf("Only aif=\"%s\" or aif=\"%s\" or aif=\"%s\" are acceptable AIFs for model=\"weinmann\" or model=\"extended\"", aif.names[1], aif.names[2], aif.names[3]), call.=FALSE)
           }
         },
         kety.orton.exp = ,
         orton.exp = {
           aif <- ifelse(is.null(aif), "orton.exp", aif)
           if (! aif %in% c("orton.exp","user")) {
             stop("Only aif=\"orton.exp\" or aif=\"user\" are acceptable AIFs for model=\"orton.exp\" or model=\"kety.orton.exp\"", call.=FALSE)
           }
         },
         kety.orton.cos= ,
         orton.cos = {
           aif <- ifelse(is.null(aif), "orton.cos", aif)
           if (! aif %in% c("orton.cos","user")) {
             stop("Only aif=\"orton.cos\" or aif=\"user\" are acceptable AIFs for model=\"orton.cos\" or model=\"kety.orton.cos\"", call.=FALSE)
           }
         },
         stop(paste("Unknown model:", model), call.=FALSE))
  p <- aifParameters(aif, user)
  if (aif == "empirical") {
    ## paste model and "empirical" together
    model <- paste(model, aif, sep=".")
  }
  func.model <- compartmentalModel(model)
  func <- function(theta, signal, time, ...) {
    out <- signal - func.model(time, theta, p)
    out[!is.na(out)]
  }
  nvoxels <- sum(img.mask)
  switch(model,
         weinmann = ,
         weinmann.empirical = ,
         kety.orton.exp = ,
         kety.orton.cos = {
           if (is.null(guess)) {
             guess <- c("th1"=-1, "th3"=-1)
           } else {
             if (length(guess) != 2 || ! all(names(guess) %in% c("th1","th3"))) {
               stop("Names of starting parameters must be \"th1\" and \"th3\"")
             }
           }
         },
         extended = ,
         extended.empirical = ,
         orton.exp = ,
         orton.cos = {
           vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
           if (is.null(guess)) {
             guess <- c("th0"=-3, "th1"=-1, "th3"=-1)
           } else {
             if (length(guess) != 3 ||
                 ! all(names(guess) %in% c("th0","th1","th3"))) {
               stop("Names of starting parameters must be \"th0\", \"th1\" and \"th3\"")
             }
           }
         },
         stop("Model/AIF combination is not supported."))

  I <- nrow(conc)
  J <- ncol(conc)
  K <- oro.nifti::nsli(conc)
  if (!is.numeric(dim(conc))) {
    I <- J <- K <- 1
  } else {
    if (length(dim(conc)) == 2)
      J <- K <- 1
  }
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  img.mask <- ifelse(img.mask > 0, TRUE, FALSE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0
  conc.list <- vector("list", nvoxels)
  for (k in 1:nvoxels) {
    conc.list[[k]] <- conc.mat[k,]
  }
  rm(conc.mat) ; gc()
  if (verbose) {
    cat("  Estimating the kinetic parameters...", fill=TRUE)
  }
  if (multicore) {
    lm.list <- parallel::mclapply(conc.list, function(x) {
      minpack.lm::nls.lm(par=guess, fn=func, control=control, signal=x, time=time, p=p)
    })
  } else {
    lm.list <- lapply(conc.list, function(x) {
      minpack.lm::nls.lm(par=guess, fn=func, control=control, signal=x, time=time, p=p)
    })
  }
  rm(conc.list) ; gc()
  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  sse <- rep(NA, nvoxels)
  for (k in 1:nvoxels) {
    if (lm.list[[k]]$info > 0 && lm.list[[k]]$info < 5) {
      ktrans$par[k] <- exp(lm.list[[k]]$par["th1"])
      kep$par[k] <- exp(lm.list[[k]]$par["th3"])
      ktrans$error[k] <- sqrt(lm.list[[k]]$hessian["th1","th1"])
      kep$error[k] <- sqrt(lm.list[[k]]$hessian["th3","th3"])
      sse[k] <- lm.list[[k]]$deviance
      if (model %in% c("extended", "orton.exp", "orton.cos",
                       "extended.empirical")) {
        vp$par[k] <- exp(lm.list[[k]]$par["th0"])
        vp$error[k] <- sqrt(lm.list[[k]]$hessian["th0","th0"])
      }
    }
  }
  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par
  B[img.mask] <- ktrans$error
  R <- list(ktrans=A, ktranserror=B, time=time)
  rm(A,B)
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par
  B[img.mask] <- kep$error
  R$kep <- A
  R$keperror <- B
  R$ve <- R$ktrans / R$kep
  rm(A,B)
  if (model %in% c("extended", "orton.exp", "orton.cos",
                   "extended.empirical")) {
    A <- B <- array(NA, c(I,J,K))
    A[img.mask] <- vp$par
    B[img.mask] <- vp$error
    R$vp <- A
    R$vperror <- B
    rm(A,B)
  }
  A <- array(NA, c(I,J,K))
  A[img.mask] <- sse
  R$sse <- A
  rm(A)
  return(R)
}
