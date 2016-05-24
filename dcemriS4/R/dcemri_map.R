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
## $Id: dcemri_map.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## setGeneric("dcemri.map")
#############################################################################
#' Pharmacokinetic Modeling of Dynamic Contrast-Enhanced MRI Data
#' 
#' Maximum-a-posteriori (MAP) estimation for single compartment models is
#' performed using literature-based or user-specified arterial input functions.
#' 
#' Implements \emph{maximum a posteriori} (MAP) estimation for the Bayesian
#' model in Schmid \emph{et al.} (2006).
#' 
#' @aliases dcemri.map
#' @param conc Matrix or array of concentration time series (last dimension
#' must be time).
#' @param time Time in minutes.
#' @param img.mask Mask matrix or array. Voxels with \code{mask=0} will be
#' excluded.
#' @param model is a character string that identifies the type of compartmental
#' model to be used.  Acceptable models include: 
#' \describe{
#' \item{\dQuote{weinmann}}{Tofts & Kermode AIF convolved with single 
#' compartment model} 
#' \item{\dQuote{extended}}{Weinmann model extended with additional vascular
#' compartment (default)} 
#' \item{\dQuote{orton.exp}}{Extended model using Orton's exponential AIF} 
#' \item{\dQuote{kety.orton.exp}}{Kety model using Orton's exponential AIF} 
#' \item{\dQuote{orton.cos}}{Extended model using Orton's raised cosine AIF} 
#' \item{\dQuote{kety.orton.cos}}{Kety model using Orton's raised cosine AIF} 
#' }
#' @param aif is a character string that identifies the parameters of the type
#' of arterial input function (AIF) used with the above model.  Acceptable
#' values are: \code{tofts.kermode} (default) or \code{fritz.hansen} for the
#' \code{weinmann} and \code{extended} models; \code{orton.exp} (default) or
#' \code{user} for the \code{orton.exp} model and \code{orton.exp} model;
#' \code{user} for the \code{orton.cos} model and \code{orton.cos} model.
#' @param user Vector of AIF parameters.  For Tofts and Kermode: \eqn{a_1},
#' \eqn{m_1}, \eqn{a_2}, \eqn{m_2}; for Orton \emph{et al.}: \eqn{A_b},
#' \eqn{\mu_b}, \eqn{A_g}, \eqn{\mu_g}.
#' @param ab.ktrans Mean and variance parameter for Gaussian prior on
#' \eqn{\log(K^{trans})}.
#' @param ab.kep Mean and variance parameter for Gaussian prior on
#' \eqn{\log(k_{ep})}.
#' @param ab.vp Hyper-prior parameters for the Beta prior on \eqn{v_p}{vp}.
#' @param ab.tauepsilon Hyper-prior parameters for observation error Gamma
#' prior.
#' @param maxit The maximum number of iterations for the optimization
#' procedure.
#' @param samples If \code{TRUE} output includes samples drawn from the
#' posterior distribution for all parameters.
#' @param multicore If \code{TRUE} algorithm is parallelized using
#' \pkg{multicore}.
#' @param verbose Logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @param ... Additional parameters to the function.
#' @return Parameter estimates and their standard errors are provided for the
#' masked region of the multidimensional array.  The multi-dimensional arrays
#' are provided in \code{nifti} format.
#' 
#' They include: 
#' \item{ktrans}{Transfer rate from plasma to the extracellular,
#' extravascular space (EES).} 
#' \item{kep}{Rate parameter for transport from the EES to plasma.} 
#' \item{ve}{Fractional occupancy by EES (the ratio between ktrans and kep).} 
#' \item{vp}{Fractional occupancy by plasma.}
#' \item{sigma2}{The residual sum-of-squares from the model fit.}
#' \item{time}{Acquisition times (for plotting purposes).
#' } 
#' Note, not all parameters are available under all models choices.
#' @author Volker Schmid \email{volkerschmid@@users.sourceforge.net}
#' @seealso \code{\link{dcemri.lm}}, \code{\link{dcemri.bayes}}
#' @references 
#' Schmid, V., Whitcher, B., Padhani, A.R., Taylor, N.J. and Yang,
#' G.-Z.  (2006) Bayesian methods for pharmacokinetic models in dynamic
#' contrast-enhanced magnetic resonance imaging, \emph{IEEE Transactions on
#' Medical Imaging}, \bold{25} (12), 1627-1636.
#' @keywords models
#' @examples
#' 
#' data("buckley")
#' xi <- seq(5, 300, by=5)
#' img <- array(t(breast$data)[,xi], c(13,1,1,60))
#' mask <- array(TRUE, dim(img)[1:3])
#' time <- buckley$time.min[xi]
#' 
#' ## MAP estimation with extended Kety model and Fritz-Hansen default AIF
#' fit.map.vp <- dcemri.map(img, time, mask, aif="fritz.hansen")
#' ## Nonlinear regression with extended Kety model and Fritz-Hansen default AIF
#' fit.lm.vp <- dcemri.lm(img, time, mask, aif="fritz.hansen")
#' 
#' plot(breast$ktrans, fit.map.vp$ktrans, xlim=c(0,1), ylim=c(0,1),
#'      xlab=expression(paste("True ", K^{trans})),
#'      ylab=expression(paste("Estimated ", K^{trans})))
#' points(breast$ktrans, fit.lm.vp$ktrans, pch=3)
#' abline(0, 1, lwd=2, col=2)
#' legend("bottomright", c("MAP Estimation (fritz.hansen)",
#'                         "Levenburg-Marquardt (fritz.hansen)"), pch=c(1,3))
#' 
#' ## MAP estimation with Kety model and Fritz-Hansen default AIF
#' fit.map <- dcemri.map(img, time, mask, model="weinmann", aif="fritz.hansen")
#' ## Nonlinear regression with Kety model and Fritz-Hansen default AIF
#' fit.lm <- dcemri.lm(img, time, mask, model="weinmann", aif="fritz.hansen")
#' 
#' cbind(breast$kep, fit.lm$kep[,,1], fit.lm.vp$kep[,,1], fit.map$kep[,,1],
#'       fit.map.vp$kep[,,1])
#' cbind(breast$ktrans, fit.lm$ktrans[,,1], fit.lm.vp$ktrans[,,1],
#'       fit.map$ktrans[,,1], fit.map.vp$ktrans[,,1])
#' 
#' @export
#' @docType methods
#' @rdname dcemri.map-methods
setGeneric("dcemri.map", function(conc, ...) standardGeneric("dcemri.map"))
#' @export
#' @rdname dcemri.map-methods
#' @aliases dcemri.map,array-method
setMethod("dcemri.map", signature(conc="array"),
          function(conc, time, img.mask, model="extended", aif=NULL,
                   user=NULL, ab.ktrans=c(0,1), ab.kep=ab.ktrans,
                   ab.vp=c(1,19), ab.tauepsilon=c(1,1/1000), maxit=5000,
                   samples=FALSE, multicore=FALSE, verbose=FALSE, ...)
          .dcemriWrapper("dcemri.map", conc, time, img.mask, model,
                         aif, user, ab.ktrans, ab.kep, ab.vp,
                         ab.tauepsilon, maxit, samples, multicore, verbose,
                         ...))

.dcemri.map.single <- function(conc, time, posterior, parameter,
                               transform, start, hyper, aif, maxit,
                               verbose=FALSE) {
  if (any(is.na(conc))) {
    return(NA)
  } else {
    map <- optim(par=start, fn=posterior, conc=conc, time=time,
                 hyper=hyper, aif=aif, control=list("maxit"=maxit))
    p <- length(parameter)
    return.list <- list()
    for (i in 1:p) {
      return.list[[parameter[i]]] <- transform[[i]](map$par[i])
    }
    return(return.list)
  }
}

.dcemri.map <- function(conc, time, img.mask, model="extended", aif=NULL,
                        user=NULL, ab.ktrans=c(0,1), ab.kep=ab.ktrans,
                        ab.vp=c(1,19), ab.tauepsilon=c(1,1/1000),
                        maxit=5000, samples=FALSE, multicore=FALSE,
                        verbose=FALSE, ...) {
  switch(model,
         weinmann = ,
         extended = {
           aif <- ifelse(is.null(aif), "tofts.kermode", aif)
           if (! aif %in% c("tofts.kermode","fritz.hansen"))
             stop("Only aif=\"tofts.kermode\" or aif=\"fritz.hansen\" are acceptable AIFs for model=\"weinmann\" or model=\"extended\"", call.=FALSE)
         },
         orton.exp = ,
         kety.orton.exp = {
           aif <- ifelse(is.null(aif), "orton.exp", aif)
           if (! aif %in% c("orton.exp","user"))
             stop("Only aif=\"orton.exp\" or aif=\"user\" are acceptable AIFs for model=\"orton.exp\" or model=\"kety.orton.exp\"", call.=FALSE)
         },
         orton.cos = ,
         kety.orton.cos = {
           aif <- ifelse(is.null(aif), "orton.cos", aif)
           if (! aif %in% c("orton.cos","user"))
             stop("Only aif=\"orton.cos\" or aif=\"user\" are acceptable AIFs for  model=\"orton.cos\" or model=\"kety.orton.cos\"", call.=FALSE)
         },
         stop(paste("Unknown model:",model), call.=FALSE))

  I <- nrow(conc)
  J <- ncol(conc)
  K <- oro.nifti::nsli(conc)

  if (!is.numeric(dim(conc))) {
    I <- J <- K <- 1
  } else {
    if (length(dim(conc)) == 2) {
      J <- K <- 1
    }
    if (length(dim(conc)) == 3) {
      K <- 1
    }
  }

  aif.parameter <- aifParameters(aif, user)
  func.model <- compartmentalModel(model)
  inverse <- function(x) {
    return(1/x)
  }
  ident <- function(x) {
    return(x)
  }

  switch(model,
         weinmann =,
         kety.orton.exp = ,
         kety.orton.cos = {
           parameter <- c("ktrans", "kep", "sigma2")
           transform <- c(exp, exp, inverse)
           hyper <- c(ab.ktrans, ab.kep, ab.tauepsilon)
           start <- c(exp(hyper[1]), exp(hyper[3]), hyper[5] * hyper[6])
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             tauepsilon <- par[3]
             conc.hat <- func.model(time, c(gamma, theta), aif)
             p <- (log(dnorm(gamma, hyper[1], hyper[2])) +
                   log(dnorm(theta, hyper[3], hyper[4])) +
                   log(dgamma(tauepsilon, hyper[5], rate=hyper[6])) +
                   sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon)))))
             p <- ifelse(is.na(p), 1e-6, p)
             return(-p)

           }
	 },
         extended = {
           inverse <- function(x) {
             1/x
           }
           ident <- function(x) {
             x
           }
           parameter <- c("ktrans", "kep", "vp", "sigma2")
           transform <- c(exp, exp, exp, inverse)
           hyper <- c(ab.ktrans, ab.kep, ab.vp, ab.tauepsilon)
           start <- c(exp(hyper[1]), exp(hyper[3]),
                      log(hyper[5]/(hyper[5]+hyper[6])), hyper[7]*hyper[8])
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             vp <- par[3]
             tauepsilon <- par[4]
             T <- length(time)
             conc.hat <- func.model(time, c(vp, gamma, theta), aif)
                         #ifelse(time > 0,
                         #       (exp(vp) * extraterm(time, aif) + exp(gamma) *
                         #        convterm(exp(theta), time, aif)),
                         #       0)
             p <- (log(dnorm(gamma, hyper[1], hyper[2])) +
                   log(dnorm(theta, hyper[3], hyper[4])) +
                   log(dgamma(tauepsilon, hyper[7], rate=hyper[8])) +
                   log(dbeta(exp(vp),hyper[5], hyper[6])) +
                   sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon)))))
             p <- ifelse(is.na(p), -1e-6, p)
             return(-p)
           }
	 },
         orton.exp = {
           inverse <- function(x) {
             1/x
           }
           ident <- function(x) {
             x
           }
           parameter <- c("ktrans", "kep", "vp", "sigma2")
           transform <- c(exp, exp, ident, inverse)
           hyper <- c(ab.ktrans, ab.kep, ab.vp, ab.tauepsilon)
           start <- c(exp(hyper[1]), exp(hyper[3]),
                      hyper[5]/(hyper[5]+hyper[6]), hyper[7]*hyper[8])
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             vp <- par[3]
             tauepsilon <- par[4]
             T <- length(time)
             conc.hat <- func.model(time, c(vp, gamma, theta), aif)
             p <- (log(dnorm(gamma, hyper[1], hyper[2])) +
                   log(dnorm(theta, hyper[3], hyper[4])) +
                   log(dgamma(tauepsilon, hyper[7], rate=hyper[8])) +
                   log(dbeta(vp, hyper[5], hyper[6])) +
                   sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon)))))
             p <- ifelse(is.na(p), 1e-6, p)
             return(-p)
           }
	 },
         kety.orton.exp = {
           inverse <- function(x) {
             1/x
           }
           parameter <- c("ktrans", "kep", "sigma2")
           transform <- c(exp, exp, inverse)
           hyper <- c(ab.ktrans, ab.kep, ab.tauepsilon)
           start <- c(exp(hyper[1]), exp(hyper[3]), hyper[5]*hyper[6])
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             tauepsilon <- par[3]
             T <- length(time)
             conc.hat <- func.model(time, c(gamma, theta), aif)
             p <- (log(dnorm(gamma, hyper[1], hyper[2])) +
                   log(dnorm(theta, hyper[3], hyper[4])) +
                   log(dgamma(tauepsilon, hyper[5], rate=hyper[6])) +
                   sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon)))))
             p <- ifelse(is.na(p), 1e-6, p)
             return(-p)
           }
	 },
        orton.cos = {
           inverse <- function(x) {
             1/x
           }
           ident <- function(x) {
             x
           }
           parameter <- c("ktrans", "kep", "vp", "sigma2")
           transform <- c(exp, exp, exp, inverse)
           hyper <- c(ab.ktrans, ab.kep, ab.vp, ab.tauepsilon)
           start <- c(hyper[1], hyper[3],
                      log(hyper[5] / (hyper[5] + hyper[6])),
                      hyper[7] * hyper[8])
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             theta0 <- par[3]
             tauepsilon <- par[4]
             T <- length(time)
             conc.hat <- func.model(time, c(theta0, gamma, theta), aif)
             p <- (log(dnorm(gamma, hyper[1], hyper[2])) +
                   log(dnorm(theta, hyper[3], hyper[4])) +
                   log(dgamma(tauepsilon, hyper[7], rate=hyper[8])) +
                   log(dbeta(exp(theta0), hyper[5], hyper[6])) +
                   sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon)))))
             p <- ifelse(is.na(p), 1e-6, p)
             return(-p)
           }
	 },
         stop("Model is not supported."))

  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  nvoxels <- sum(img.mask)
  img.mask <- ifelse(img.mask > 0, TRUE, FALSE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0
  conc.list <- vector("list", nvoxels) # list()
  for (i in 1:nvoxels) {
    conc.list[[i]] <- conc.mat[i,]
  }

  if (verbose) {
    cat("  Estimating the kinetic parameters...", fill=TRUE)
  }
  if (multicore) {
    fit <- parallel::mclapply(conc.list, FUN=.dcemri.map.single, time=time,
                    posterior=posterior, parameter=parameter,
                    transform=transform, start=start, hyper=hyper,
                    aif=aif.parameter, maxit=maxit, verbose=verbose)
  } else {
    fit <- lapply(conc.list, FUN=.dcemri.map.single, time=time,
                  posterior=posterior, parameter=parameter,
		  transform=transform, start=start, hyper=hyper,
                  aif=aif.parameter, maxit=maxit, verbose=verbose)
  }

  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  sigma2 <- rep(NA, nvoxels)
  vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  for (k in 1:nvoxels) {
    try(ktrans$par[k] <- fit[[k]]$ktrans, silent=TRUE)
    try(kep$par[k] <- fit[[k]]$kep, silent=TRUE)
    if (model %in% c("extended", "orton.exp", "orton.cos")) {
      try(vp$par[k] <- fit[[k]]$vp, silent=TRUE)
    }
    try(sigma2[k] <- fit[[k]]$sigma2, silent=TRUE)
  }
  A <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par

  returnable <- list(ktrans=A)
  A <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par

  returnable[["kep"]] <- A
  if (model %in% c("extended", "orton.exp", "orton.cos")) {
    A <- array(NA, c(I,J,K))
    A[img.mask] <- vp$par
    returnable[["vp"]] <- A
  }
  A <- array(NA, c(I,J,K))
  A[img.mask] <- sigma2

  returnable[["sigma2"]] <- A
  returnable[["time"]] <- time
  returnable[["ve"]] <- returnable$ktrans / returnable$kep

  return(returnable)
}
