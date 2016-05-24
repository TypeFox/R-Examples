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
## $Id: dcemri_spline.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## setGeneric("dcemri.spline")
#############################################################################
#' Bayesian P-Splines for Dynamic Contrast-Enhanced MRI Data
#' 
#' Quantitative analysis of DCE-MRI typically involves the convolution of an
#' arterial input function (AIF) with a nonlinear pharmacokinetic model of the
#' contrast agent concentration.  This function takes a semi-parametric
#' penalized spline smoothing approach, with which the AIF is convolved with a
#' set of B-splines to produce a design matrix using locally adaptive smoothing
#' parameters based on Bayesian penalized spline models (P-splines).
#' 
#' See Schmid \emph{et al.} (2009) for more details.
#' 
#' @aliases dcemri.spline dcemri.spline,array-method dcemri.spline.single
#' @param conc Matrix or array of concentration time series (last dimension
#' must be time).
#' @param ... Additional variables defined by the method.  
#' @param time Time in minutes.
#' @param img.mask Mask matrix or array. Voxels with \code{mask = 0} will be
#' excluded.
#' @param time.input Time in minutes for observed arterial input function
#' (default = \sQuote{time}).
#' @param aif is a character string that identifies the parameters of the
#' arterial input function.  Acceptable values are: \code{tofts.kermode},
#' \code{fritz.hansen} or \code{observed}.  If \code{observed} you must provide
#' the observed concentrations in \code{aif.observed}.
#' @param user \ldots{}
#' @param aif.observed is the user-defined vector of arterial concentrations
#' observed at \code{time.input} (only for \sQuote{aif}=observed).
#' @param multicore (logical) use the \pkg{parallel} package.
#' @param verbose (logical) allows text-based feedback during execution of the
#' function (default = \code{FALSE}).
#' @param samples If \code{TRUE} output includes samples drawn from the
#' posterior distribution for all parameters.
#' @param nlr If \code{TRUE}, a response model is fitted to the estimated
#' response function.
#' @param model Only if \code{nlr = TRUE} Response model fitted to the
#' estimated response function.  Acceptable values include: \code{"AATH"} or
#' \code{"weinmann"} (default).
#' @param ab.hyper Hyper priors for adaptive smoothness parameter
#' @param ab.tauepsilon Hyper-prior parameters for observation error Gamma
#' prior.
#' @param p Number of knots of B-Spline basis.
#' @param t0.compute If \code{TRUE}, the onset time will be estimated from
#' response function.
#' @param k Order of B-Splines.
#' @param knots Vector of knots.  Use this if you need unequally spaced knots.
#' @param rw Order of random walk prior.  Acceptable values are 1 and 2.
#' @param nriters Total number of iterations.
#' @param thin Thining factor.
#' @param burnin Number of iterations for burn-in.
#' @param response If \code{TRUE}, the response functions per voxel are
#' returned.
#' @param fitted If \code{TRUE}, then fitted time curved per voxel are
#' returned.
#' @return The maximum of the response function \code{Fp} for the masked region
#' is provided by default.  Where appropriate, response functions, fitted
#' functions, and parameter estimates (along with their standard errors) are
#' provided. All multi-dimensional arrays are provided in \code{nifti} format.
#' @author Volker Schmid \email{volkerschmid@@users.sourceforge.net}
#' @seealso \code{\link{dcemri.bayes}}, \code{\link{dcemri.lm}},
#' \code{\link{dcemri.map}}
#' @references 
#' Schmid, V., Whitcher, B., Padhani, A.R. and G.-Z. Yang (2009) A
#' semi-parametric technique for the quantitative analysis of dynamic
#' contrast-enhanced MR images based on Bayesian P-splines, \emph{IEEE
#' Transactions on Medical Imaging}, \bold{28} (6), 789-798.
#' @keywords models
#' @examples
#' 
#' data("buckley")
#' xi <- seq(5, 300, by=5)
#' img <- array(t(breast$data)[,xi], c(13,1,1,60))
#' mask <- array(TRUE, dim(img)[1:3])
#' time <- buckley$time.min[xi]
#' 
#' ## Generate AIF params using the orton.exp function from Buckley's AIF
#' aif <- buckley$input[xi]
#' 
#' fit.spline <- dcemri.spline(img, time, mask, aif="fritz.hansen",
#'                             nriters=125, thin=3, burnin=25, nlr=TRUE)
#' fit.spline.aif <- dcemri.spline(img, time, mask, aif="observed",
#'                                 aif.observed=aif, nriters=125, thin=3,
#'                                 burnin=25, nlr=TRUE)
#' 
#' plot(breast$ktrans, fit.spline$ktrans, xlim=c(0,1), ylim=c(0,1),
#'      xlab=expression(paste("True ", K^{trans})),
#'      ylab=expression(paste("Estimated ", K^{trans})))
#' points(breast$ktrans, fit.spline.aif$ktrans, pch=2)
#' abline(0, 1, lwd=1.5, col="red")
#' legend("right", c("fritz.hansen", "observed"), pch=1:2)
#' 
#' @export
#' @docType methods
#' @rdname dcemri.spline-methods
setGeneric("dcemri.spline", function(conc, ...) standardGeneric("dcemri.spline"))
#' @export
#' @rdname dcemri.spline-methods
#' @aliases dcemri.spline,array-method
#' @useDynLib dcemriS4 dce_spline_run
setMethod("dcemri.spline", signature(conc="array"),
          function(conc, time, img.mask, time.input=time,
                   model="weinmann", aif="tofts.kermode",
                   user=NULL, aif.observed=NULL, nriters=500,
                   thin=5, burnin=100, ab.hyper=c(1e-5,1e-5),
                   ab.tauepsilon=c(1,1/1000), k=4, p=25, rw=2,
                   knots=NULL, nlr=FALSE, t0.compute=FALSE,
                   samples=FALSE, multicore=FALSE, verbose=FALSE,
		   response=FALSE, fitted=FALSE, ...)
          .dcemriWrapper("dcemri.spline", conc, time, img.mask, time.input,
                         model, aif, user, aif.observed, nriters, thin,
                         burnin, ab.hyper, ab.tauepsilon, k, p, rw, knots,
                         nlr, t0.compute, samples, multicore, verbose,
                         response, fitted, ...))

.dcemri.spline.single <- function(conc, time, D, time.input, p, rw, knots,
                                  k, A, t0.compute=FALSE, nlr=FALSE,
                                  nriters=500, thin=5, burnin=100,
                                  ab.hyper=c(1e-5,1e-5),
                                  ab.tauepsilon=c(1,1/1000), silent=0,
                                  multicore=FALSE, model=NULL,
                                  model.func=NULL, model.guess=NULL,
                                  samples=FALSE, B=NULL) {

  ## Sanity check: conc must not contain any missing values
  if (any(is.na(conc))) {
    stop("Concentration time curves must not contain missing values.")
  }

  samplesize <- floor(nriters/thin)

  T <- length(time)
  TT <- T*T
  tau <- matrix(1000, p-rw, samplesize) # array(1000, c(p-rw,samplesize))
  beta <- matrix(0, p, samplesize) # array(0, c(p,samplesize))
  ##MAX <- array(0, c(samplesize))
  tauepsilon <- rep(1000, samplesize) # array(1000, c(samplesize))
  burnin <- min(burnin, samplesize)

  result <- .C("dce_spline_run",
               as.integer(1),
               as.integer(burnin),
               as.integer(c(1,1,1,T)),
               as.double(conc),
               as.double(tau),
               as.double(tauepsilon),
               as.double(D),
               as.integer(rw),
               as.double(beta),
               as.double(c(ab.hyper[1], ab.hyper[2], ab.tauepsilon[1],
                           ab.tauepsilon[2])),
               as.integer(p),
               as.double(1:T),
               as.double(1:T),
               as.double(1:T),
               as.double(1:TT),
               as.double(1:TT),
               as.double(1:TT),
               as.double(1:TT),
               as.double(t(D)),
               as.integer(silent),
               PACKAGE="dcemriS4")

  tau <- matrix(result[[5]], p-rw, samplesize)
  tauepsilon <- result[[6]]
  beta <- matrix(result[[9]], p, samplesize)

  ##
  ## Why is this run twice?!?
  ##

  result <- .C("dce_spline_run",
               as.integer(samplesize),
               as.integer(thin),
               as.integer(c(1,1,1,T)),
               as.double(conc),
               as.double(tau),
               as.double(tauepsilon),
               as.double(D),
               as.integer(rw),
               as.double(beta),
               as.double(c(ab.hyper[1], ab.hyper[2], ab.tauepsilon[1],
                           ab.tauepsilon[2])),
               as.integer(p),
               as.double(1:T),
               as.double(1:T),
               as.double(1:T),
               as.double(1:TT),
               as.double(1:TT),
               as.double(1:TT),
               as.double(1:TT),
               as.double(t(D)),
               as.integer(silent),
               PACKAGE="dcemriS4")

  tau <- matrix(result[[5]], p-rw, samplesize)
  tauepsilon <- result[[6]]
  beta <- matrix(result[[9]], p, samplesize)

  t0 <- 0
  if (t0.compute) {
    d <- matrix(NA, T, samplesize)
    for (j in 1:samplesize) {
      d[,j] <- D %*% beta[,j]
    }
    d1 <- apply(d, 1, quantile, probs=c(0.005), na.rm=TRUE)
    d2 <- apply(d, 1, median, na.rm=TRUE)
    du <- min(which(d1 > 0))
    ## beta.abl <- beta.abl2 <- rep(0, p)
    B2 <- splines::splineDesign(knots, time.input, k-2)
    B2 <- B2[,(1:p)+1]

    for (j in 1:samplesize) {
      beta.abl <- 1:p
      for (q in 1:(p - 1)) {
	beta.abl[q] <- (beta[q+1, j] - beta[q, j]) *
          k / (knots[q+k+1] - knots[q+1])
      }
      beta.abl[p] <- 0
      ABL2 <- A %*% B2 %*% beta.abl
      du2 <- time[du] - d2[du] / ABL2[du]
      t0[j] <- du2
    }
    if (sum(!is.na(t0)) == 0) {
      t0 <- 0
    } else {
      t0 <- median(t0, na.rm=TRUE)
    }
    if (t0 < min(time)) {
      t0 <- 0
    }
    if (t0 > max(time)) {
      t0 <- 0
    }
  }

  fitted <- vector("list", samplesize) # list()
  for (i in 1:samplesize) {
    fitted[[i]] <- B %*% beta[,i]
  }

  ##if (multicore && require("parallel")) {
  ##  MAX <- unlist(mclapply(fitted, max))
  ##} else {
    MAX <- unlist(lapply(fitted, max))
  ##}
  ##MAX <- rep(NA, samplesize)
  ##for (i in 1:samplesize) {
  ##  MAX[i] <- MAX0[[i]]
  ##}
  parameters <- list()
  if (nlr) {
    if (model=="AATH") {
      model.guess[2] <- median(MAX, na.rm=TRUE)
    }
    fcn <- function(p, time, x, N.Err, fcall, jcall)
      (x - do.call("fcall", c(list(time=time), as.list(p))))

    nls.lm.single <- function(fitted, par, fn, fcall, model, time) {
      fcall2 <- fcall
      if (length(fcall) > 1) {
	fcall2 <- fcall[[1]]
      }
      fit <- minpack.lm::nls.lm(par=par, fn=fn, fcall=fcall2, time=time, x=fitted,
                    N.Err=sqrt(300), control=list(nprint=0, ftol=10^-20))
      if (model=="AATH" && fit$par$TC < 1e-6) {
	fit <- minpack.lm::nls.lm(par=par[-3], fn=fn, fcall=fcall[[2]], time=time,
                      x=fitted, N.Err=sqrt(300),
                      control=list(nprint=0, ftol=10^-20))
	fit$par$TC <- 0
      }
      return(fit)
    }

    if (multicore) {
      out <- parallel::mclapply(fitted, nls.lm.single, par=model.guess,
                                fn=fcn, fcall=model.func, model=model,
                                time=time-t0)
    } else {
      out <- lapply(fitted, nls.lm.single, par=model.guess,
                         fn=fcn, fcall=model.func, model=model,
                         time=time-t0)
    }

    if (model=="AATH") {
      E <- F <- TC <- ve <- rep(NA, samplesize) # NULL
      for (i in 1:samplesize) {
	E[i] <- exp(out[[i]]$par$logE)
	F[i] <- exp(out[[i]]$par$logF)
	TC[i] <- out[[i]]$par$TC
	ve[i] <- exp(out[[i]]$par$logve)
      }
      parameters <- list("E"=median(E, na.rm=TRUE),
                         "F"=median(F, na.rm=TRUE),
                         "TC"=median(TC, na.rm=TRUE),
                         "ve"=median(ve, na.rm=TRUE))
      if (samples) {
	parameters$E.samples <- E
        parameters$F.samples <- F
        parameters$TC.samples <- TC
        parameters$ve.samples <- ve
      }
    }
    if (model == "weinmann") {
      ktrans <- kep <- rep(NA, samplesize) # NULL
      for (i in 1:samplesize) {
	ktrans[i] <- exp(out[[i]]$par$logktrans)
	kep[i] <- exp(out[[i]]$par$logkep)
      }
      parameters <- list(ktrans=median(ktrans, na.rm=TRUE),
                         kep=median(kep, na.rm=TRUE))
      if (samples) {
	parameters$ktrans.sample <- ktrans
        parameters$kep.sample <- kep
      }
    }
  }

  list(beta=beta, tau=tau, tauepsilon=tauepsilon, t0=t0,
       Fp=median(MAX, na.rm=TRUE), Fp.samples=MAX, fitted=fitted,
       par=parameters)
}

.dcemri.spline <- function(conc, time, img.mask, time.input=time,
                          model="weinmann", aif="tofts.kermode",
                          user=NULL, aif.observed=NULL, nriters=500,
                          thin=5, burnin=100, ab.hyper=c(1e-5,1e-5),
                          ab.tauepsilon=c(1,1/1000), k=4, p=25, rw=2,
                          knots=NULL, nlr=FALSE, t0.compute=FALSE,
                          samples=FALSE, multicore=FALSE, verbose=FALSE,
                          response=FALSE, fitted=FALSE, ...) {

  ## dcemri.spline - a function for fitting Bayesian Penalty Splines to
  ## DCE-MRI images and computing kinetic parameters
  ##
  ## authors: Volker Schmid, Brandon Whitcher
  ##
  ## input:
  ##        conc: array of Gd concentration,
  ##        time: timepoints of aquisition,
  ##        img.mask: array of voxels to fit,
  ##        D(=0.1): Gd dose in mmol/kg,
  ##        model: AIF... "weinman" or "parker",
  ##
  ## output: list with ktrans, kep, ve, std.error of ktrans and kep
  ##         (ktranserror and keperror)
  ##

  ## function to make precision matrix for random walk
  R <- function(taux, rw) {
    RR <- matrix(0, length(taux)+rw, length(taux)+rw)
    if (rw == 0) {
      for (i in 1:length(taux)) {
	RR[i,i] <- taux[i]
      }
    }
    if (rw == 1) {
      for (i in 1:length(taux)) {
	RR[i,i] <- RR[i,i] + taux[i]
	RR[i+1,i+1] <- RR[i+1,i+1] + taux[i]
	RR[i+1,i] <- RR[i+1,i] - taux[i]
	RR[i,i+1] <- RR[i,i+1] - taux[i]
      }
    }
    if (rw == 2) {
      for (i in 1:length(taux)) {
	RR[i,i] <- RR[i,i] + taux[i]
	RR[i+1,i+1] <- RR[i+1,i+1] + 4*taux[i]
	RR[i+2,i+2] <- RR[i+2,i+2] + taux[i]
	RR[i+1,i] <- RR[i+1,i] - 2*taux[i]
	RR[i,i+1] <- RR[i,i+1] - 2*taux[i]
	RR[i+2,i+1] <- RR[i+2,i+1] - 2*taux[i]
	RR[i+1,i+2] <- RR[i+1,i+2] - 2*taux[i]
	RR[i+2,i] <- RR[i+2,i] + taux[i]
	RR[i,i+2] <- RR [i,i+2] + taux[i]
      }
    }
    return(RR)
  }

  ## Main function

  knotpoints <- p
  if (is.null(knots)) {
    a <- min(time) - 1e-3 - (k-1) * (max(time)-min(time)) / (knotpoints-k+1)
    b <- 1e-3 + max(time) + k * (max(time)-min(time)) / (knotpoints-k+1)
    knots <- seq(a, b, (2e-3 + max(time) - min(time)) / (knotpoints - k + 1))
  }
  mod <- model
  nvoxels <- sum(img.mask)
  I <- nrow(conc)
  J <- ncol(conc)
  K <- oro.nifti::nsli(conc)
  T <- length(time)
  Ti <- length(time.input)

  if (!is.numeric(dim(conc))) {
    I <- J <- K <- 1
  } else {
    if (length(dim(conc)) == 2) {
      J <- K <- 1
    } else {
      if (length(dim(conc)) == 3) {
        K <- 1
      }
    }
  }

  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  img.mask <- ifelse(img.mask > 0, TRUE, FALSE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0
  conc.list <- vector("list", nvoxels) # list()
  for (i in 1:nvoxels) {
    conc.list[[i]] <- conc.mat[i,]
  }

  switch(aif,
    tofts.kermode = {
      D <- 0.1; a1 <- 3.99; a2 <- 4.78; m1 <- 0.144; m2 <- 0.0111
      input <- D*(a1*exp(-m1*time)+a2*exp(-m2*time))
    },
    fritz.hansen = {
      D <- 1; a1 <- 2.4; a2 <- 0.62; m1 <- 3.0; m2 <- 0.016
      input <- D*(a1*exp(-m1*time)+a2*exp(-m2*time))
    },
    observed = {
      input <- aif.observed
    },
    print("WARNING: AIF parameters must be specified!"))

  model.func <- model.guess <- NULL

  if (nlr && model != "weinmann" && model != "AATH") {
    stop("Only model=\"weinmann\" or model=\"AATH\" acceptable AIFs for nlr=TRUE")
  }

  if (model == "weinmann") {
    ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
    sigma2 <- rep(NA, nvoxels)
    model.func <- function(time, logktrans, logkep) {
      ktrans <- exp(logktrans)
      kep <- exp(logkep)
      erg <- ktrans * exp(-kep*(time))
      eval(erg)
    }
    model.guess <- list(logktrans=-1, logkep=0)
  }

  if (model=="AATH") {
    E <- F <- TC <- ve <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
    sigma2 <- rep(NA, nvoxels)

    model.func <- list()
    model.func[[1]] <- function(time, logF, logE, TC, logve) {
      TC2 <- 2*TC
      if (TC < 0) {
        TC2 <- 0
      }
      E <- exp(logE)
      F <- exp(logF)
      ve <- exp(logve)
      kep <- E * F/ve
      erg <- E * exp(-kep*(time-TC))
      erg[time<TC] <- 1 # - time[time<TC2]*(1-E) / TC
      erg <- erg*F
      if (TC < 0) {
	erg <- rep(-10^16, T)
      }
      eval(erg)
    }
    model.func[[2]] <- function(time, logF, logE, logve) {
      E <- exp(logE)
      F <- exp(logF)
      ve <- exp(logve)
      kep <- E * F / ve
      erg <- E * exp(-kep * time)
      erg <- erg * F
      eval(erg)
    }
    model.guess <- list(logE=log(0.6), logF=log(2), TC=0, logve=log(0.05))
  }

  ## define A and B

  p <- length(knots) - k
  B <- splines::splineDesign(knots, time.input, k, outer.ok=TRUE)
  if (sum(B[, dim(B)[2]] == 0) == dim(B)[1]) {
    B <- B[,-dim(B)[2]]
  }
  if (sum(B[,1] == 0) == dim(B)[1]) {
    B <- B[,-1]
  }
  p <- dim(B)[2]
  A <- matrix(0, T, Ti)
  ni <- time
  for (i in 1:T) {
    for (j in 1:Ti) {
      if (time.input[j] <= time[i]) {
	ni[i] <- j
      }
    }
  }
  for (i in 1:T) {
    for (j in 1:ni[i]) {
      A[i,j] <- input[1+ni[i]-j]
    }
  }
  A <- A * mean(diff(time.input), na.rm=TRUE)
  A[is.na(A)] <- 0
  D <- A %*% B
  ## T <- length(time)

  if (verbose) {
    cat("  Estimating the parameters...", fill=TRUE)
  }

  if (multicore) {
    fit <- parallel::mclapply(conc.list, FUN=.dcemri.spline.single, time=time,
                              D=D, time.input=time.input, p=p, rw=rw, knots=knots,
                              k=k, A=A, nriters=nriters, thin=thin, burnin=burnin,
                              ab.hyper=ab.hyper, ab.tauepsilon=ab.tauepsilon,
                              t0.compute=t0.compute, nlr=nlr, multicore=TRUE,
                              model=model, model.func=model.func,
                              model.guess=model.guess, samples=samples, B=B)
  } else {
    fit <- lapply(conc.list, FUN=.dcemri.spline.single, time=time, D=D,
                  time.input=time.input, p=p, rw=rw, knots=knots, k=k,
                  A=A, nriters=nriters, thin=thin, burnin=burnin,
                  ab.hyper=ab.hyper, ab.tauepsilon=ab.tauepsilon,
                  t0.compute=t0.compute, nlr=nlr, multicore=FALSE,
                  model=model, model.func=model.func,
                  model.guess=model.guess, samples=samples, B=B)
  }

  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }

  samplesize <- floor(nriters/thin)

  t0 <- numeric(nvoxels)
  for (k in 1:nvoxels) {
    t0 <- fit[[k]]$t0
  }
  t0.img <- array(NA, c(I,J,K))
  t0.img[img.mask] <- t0
  t0 <- t0.img
  rm(t0.img)

  Fp <- numeric(nvoxels)
  for (k in 1:nvoxels) {
    Fp[k] <- fit[[k]]$Fp
  }
  Fp.img <- array(NA, c(I,J,K))
  Fp.img[img.mask] <- Fp
  Fp.samples <- matrix(NA, nvoxels, samplesize)
  for (i in 1:nvoxels) {
    Fp.samples[i,] <- fit[[k]]$Fp.samples
  }
  if (samples) {
    Fp <- array(NA, c(I,J,K,samplesize))
    for (j in 1:samplesize) {
      Fp[,,,j][img.mask] <- Fp.samples[,j]
    }
  }

  if (nlr) {
    ktrans <- ve <- rep(NA, nvoxels) # NULL
    if (model=="weinmann") {
      kep <- NULL
      if (samples) {
	ktrans.samples <- matrix(NA, nvoxels, samplesize)
      }
      for (k in 1:nvoxels) {
	ktrans[k] <- fit[[k]]$par$ktrans
      	if (samples) {
	  ktrans.samples[k,] <- fit[[k]]$par$ktrans.samples
        }
      }
      if (samples) {
	kep.samples <- array(NA, c(nvoxels,samplesize))
      }
      for (k in 1:nvoxels) {
	kep[k] <- fit[[k]]$par$kep
	if (samples) {
	  kep.samples[k,] <- fit[[k]]$par$kep.samples
        }
      }
      ve <- ktrans / kep
      if (samples) {
	ve.samples <- ktrans.samples / kep.samples
      }
    }
    if (model == "AATH") {
      E <- F <- TC <- rep(NA, nvoxels)
      if (samples) {
	E.samples <- matrix(NA, nvoxels, samplesize)
      }
      for (k in 1:nvoxels) {
	E[k] <- fit[[k]]$par$E
	if (samples) {
	  E.samples[k,] <- fit[[k]]$par$E.samples
        }
      }
      if (samples) {
	F.sample <- matrix(NA, nvoxels, samplesize)
      }
      for (k in 1:nvoxels) {
	F[k] <- fit[[k]]$par$F
	if (samples) {
	  F.samples[k,] <- fit[[k]]$par$F.samples
        }
      }
      if (samples) {
	TC.samples <- matrix(NA, nvoxels, samplesize)
      }
      for (k in 1:nvoxels) {
	TC[k] <- fit[[k]]$par$TC
	if (samples) {
	  TC.samples[k,] <- fit[[k]]$par$TC.samples
        }
      }
      if (samples) {
	ve.samples <- matrix(NA, nvoxels, samplesize)
      }
      for (k in 1:nvoxels) {
	ve[k] <- fit[[k]]$par$ve
	if (samples) {
	  ve[k,] <- fit[[k]]$par$ve.samples
        }
      }
      ktrans <- E*F
      if (samples) {
	ktrans.samples <- E.samples / F.samples
      }
    }
  }

  beta.sample <- array(NA, c(nvoxels,p,samplesize))
  response.sample <- array(NA, c(nvoxels,Ti,samplesize))
  fitted.sample <- array(NA, c(nvoxels,T,samplesize))
  for (k in 1:nvoxels) {
    beta.sample[k,,] <- fit[[k]]$beta
    for (j in 1:samplesize) {
      response.sample[k,,j] <- fit[[k]]$fitted[[j]]
      fitted.sample[k,,j] <- D %*% beta.sample[k,,j]
    }
  }

  if ((response || fitted) && samples) {
    beta <- array(NA, c(I,J,K,p,samplesize))     # a 5-dimensional array!!!
  }
  if (response || fitted) {
    beta.med <- array(NA, c(I,J,K,p))
  }
  if (fitted && samples) {
    fitted <- array(NA, c(I,J,K,T,samplesize))   # a 5-dimensional array!!!
  }
  if (fitted) {
    fitted.med <- array(NA, c(I,J,K,T))
  }
  if (response && samples) {
    response <- array(NA, c(I,J,K,T,samplesize)) # a 5-dimensional array!!!
  }
  if (response) {
    response.med <- array(NA, c(I,J,K,T))
  }
  ktrans.med <- ve.med <- array(NA, c(I,J,K))
  ktrans.med[img.mask] <- ktrans
  ve.med[img.mask] <- ve
  if (model=="weinmann") {
    kep.med <- array(NA, c(I,J,K))
    kep.med[img.mask] <- kep
  }
  if (model=="AATH") {
    E.med <- F.med <- TC.med <- array(NA, c(I,J,K))
    E.med[img.mask] <- E
    F.med[img.mask] <- F
    TC.med[img.mask] <- TC
  }
  if (samples) {
    ktrans <- ve <- array(NA, c(I,J,K,samplesize))
    for (i in 1:samplesize) {
      ktrans[,,,i][img.mask] <- ktrans.samples[,i]
      ve[,,,i][img.mask] <- ve.samples[,i]
    }
    if (model=="weinmann") {
      kep <- array(NA,c(I,J,K,samplesize))
      for (i in 1:samplesize) {
        kep[,,,i][img.mask] <- kep.samples[,i]
      }
    }
    if (model=="AATH") {
      E <- F <- TC <- array(NA, c(I,J,K,samplesize))
      for (i in 1:samplesize) {
        E[,,,i][img.mask] <- E.samples[,i]
        F[,,,i][img.mask] <- F.samples[,i]
        TC[,,,i][img.mask] <- TC.samples[,i]
      }
    }
  }
# }
  if (fitted || response) {
    if (I > 1) {
      for (j in 1:p) {
        beta.med[,,,j][img.mask] <- apply(beta.sample[,j,], 1, median)
        if (samples)
          for (i in 1:samplesize) {
            beta[,,,j,i][img.mask] <- beta.sample[,j,i]
          }
      }
    } else {
      for (j in 1:p) {
        beta.med[,,,j][img.mask] <- median(beta.sample[,j,])
        for (i in 1:samplesize) {
          beta[,,,j,i][img.mask] <- beta.sample[,j,i]
        }
      }
    }
  }

if (fitted || response) {
    if (I == 1)
    for (j in 1:T) {
      if (response)
        response.med[,,,j][img.mask] <- median(response.sample[,j,])
      if (fitted)
        fitted.med[,,,j][img.mask] <- median(fitted.sample[,j,])
      }
    if (I > 1)
    for (j in 1:T) {
      if (response)
        response.med[,,,j][img.mask] <- apply(response.sample[,j,], 1, median)
      if (fitted)
        fitted.med[,,,j][img.mask] <- apply(fitted.sample[,j,], 1, median)
      }
    if (samples)
    for (i in 1:samplesize)
    for (j in 1:T) {
      if (response)
        response[,,,j,i][img.mask] <- response.sample[,j,i]
      if (fitted)
        fitted[,,,j,i][img.mask] <- fitted.sample[,j,i]
      }
  }

  return.list <- vector("list")
  if (response || fitted) {
    return.list[["beta"]] <- beta.med
  }
  if ((response || fitted) && samples) {
    return.list[["beta.sample"]] <- beta
  }
  if (fitted) {
    return.list[["fit"]] <- fitted.med
  }
  if (fitted && samples) {
    return.list[["fit.sample"]] <- fitted
  }
  if (response) {
    return.list[["response"]] <- response.med
  }
  if (response && samples) {
    return.list[["response.sample"]] <- response
  }
  return.list[["Fp"]] <- Fp.img
  return.list[["A"]] <- A
  return.list[["B"]] <- B
  return.list[["D"]] <- D
  if (t0.compute) {
    return.list[["t0"]] <- t0
  }
  if (nlr) {
    if (model=="weinmann") {
      return.list[["kep"]] <- kep.med
    }
    if (model=="AATH") {
      return.list[["E"]] <- E.med
      return.list[["F"]] <- F.med
      return.list[["TC"]] <- TC.med
    }
    return.list[["ktrans"]] <- ktrans.med
    return.list[["ve"]] <- ve.med
    if (samples) {
      return.list[["ktrans.sample"]] <- ktrans
      return.list[["ve.sample"]] <- ve
      if (model=="weinmann") {
        return.list[["kep.samples"]] <- kep
      }
      if (model=="AATH") {
        return.list[["E.samples"]] <- E
        return.list[["F.samples"]] <- F
        return.list[["TC.samples"]] <- TC
      }
    }
  }
  if (samples) {
    return.list[["Fp.samples"]] <- Fp
  }

  return(return.list)
}
