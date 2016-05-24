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
## setGeneric("dcemri.bayes")
#############################################################################

#' Bayesian Methods for Pharmacokinetic Modeling of Dynamic Contrast-Enhanced
#' MRI Data
#' 
#' Bayesian analysis of contrast agent concentration time curves from DCE-MRI.
#' 
#' See Schmid \emph{et al.} (2006) for details.
#' 
#' @aliases dcemri.bayes dcemri.bayes,array-method dcemri.bayes.single
#' @param conc Matrix or array of concentration time series (last dimension
#' must be time).
#' @param time Time in minutes.
#' @param img.mask Mask matrix or array. Voxels with mask=0 will be excluded.
#' @param model is a character string that identifies the type of compartmental
#' model to be used.  Acceptable models include: 
#' \describe{
#' \item{"weinmann"}{Tofts & Kermode AIF convolved with single 
#' compartment model} 
#' \item{"extended"}{Weinmann model extended with additional 
#' vascular compartment (default)}
#' \item{"orton.exp"}{Extended model using Orton's exponential AIF} 
#' \item{"kety.orton.exp"}{Kety model using Orton's exponential AIF}
#' }
#' @param aif is a character string that identifies the parameters of the type
#' of arterial input function (AIF) used with the above model.  Acceptable
#' values are: 
#' \code{tofts.kermode} (default) or \code{fritz.hansen} for the
#' \code{weinmann} and \code{extended} models; 
#' \code{orton.exp} (default) or \code{user} for the \code{orton.exp} and 
#' \code{kety.orton.exp} model.
#' @param user Vector of AIF parameters.  For Tofts and Kermode: \eqn{a_1}{a1},
#' \eqn{m_1}{m1}, \eqn{a_2}{a2}, \eqn{m_2}{m2}; for Orton \emph{et al.}:
#' \eqn{A_b}{Ab}, \eqn{\mu_b}{mub}, \eqn{A_g}{Ag}, \eqn{\mu_g}{mug}.
#' @param nriters Total number of iterations.
#' @param thin Thining factor.
#' @param burnin Number of iterations for burn-in.
#' @param tune Number for iterations for tuning.  The algorithm will be tuned
#' to an acceptance rate between 0.3 and 0.6.
#' @param ab.ktrans Mean and variance parameter for Gaussian prior on
#' \eqn{\log(K^{trans})}{log(Ktrans)}.
#' @param ab.kep Mean and variance parameter for Gaussian prior on
#' \eqn{\log(k_{ep})}{log(kep)}.
#' @param ab.vp Hyper-prior parameters for the Beta prior on \eqn{v_p}{vp}.
#' @param ab.tauepsilon Hyper-prior parameters for observation error Gamma
#' prior.
#' @param samples If \code{TRUE} output includes samples drawn from the
#' posterior distribution for all parameters.
#' @param multicore If \code{TRUE} algorithm is parallelized using
#' \pkg{multicore}.
#' @param verbose Logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @param dic If \code{TRUE}, the deviance information criterion (DIC) and
#' effective number of parameters (pD) will be computed.  If 
#' \code{"samples = TRUE"}, then samples of the DIC and pD will be given.
#' @param vp Fractional occupancy in the plasma space.
#' @param ... Additional parameters to the function.
#' @return Parameter estimates and their standard errors are provided for the
#' masked region of the multidimensional array.  All multi-dimensional arrays
#' are output in \code{nifti} format.
#' They include: 
#' \item{ktrans}{Transfer rate from plasma to the extracellular,
#' extravascular space (EES).} 
#' \item{ktranserror}{Error on \eqn{K^{trans}}{ktrans}.} 
#' \item{kep}{Rate parameter for transport from the EES to plasma.} 
#' \item{keperror}{Error on \eqn{k_{ep}}{kep}.}
#' \item{ve}{Fractional occupancy by EES (the ratio between ktrans and kep).}
#' \item{vperror}{Error on \eqn{v_e}{ve}.} 
#' \item{vp}{Fractional occupancy by plasma.} 
#' \item{sigma2}{The residual sum-of-squares from the model fit.}
#' \item{time}{Acquisition times (for plotting purposes).} 
#' \item{DIC}{Deviance information criterion.} 
#' \item{DIC.map}{Contribution to DIC per voxel.}
#' \item{pD}{Effective number of parameters.} 
#' \item{pD.map}{Constribution to pD per voxel.} 
#' Note, not all parameters are available under all models choices.
#' @author Volker Schmid \email{volkerschmid@@users.sourceforge.net}
#' @seealso \code{\link{dcemri.lm}}, \code{\link{dcemri.map}},
#' \code{\link{dcemri.spline}}
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
#' ## Bayesian estimation with Fritz-Hansen default AIF
#' fit.bayes <- dcemri.bayes(img, time, mask, aif="fritz.hansen",
#'                          nriters=1000, thin=2, burnin=200)
#' 
#' ## Bayesian estimation with "orton.exp" function fit to Buckley's AIF
#' aif <- buckley$input[xi]
#' aifparams <- orton.exp.lm(time, aif)
#' aifparams$D <- 1
#' fit.bayes.aif <- dcemri.bayes(img, time, mask, model="orton.exp",
#'                               aif="user", user=aifparams,
#'                               nriters=1000, thin=2, burnin=200)
#' 
#' plot(breast$ktrans, fit.bayes$ktrans, xlim=c(0,1), ylim=c(0,1),
#'      xlab=expression(paste("True ", K^{trans})),
#'      ylab=expression(paste("Estimated ", K^{trans}, " (Bayesian)")))
#' points(breast$ktrans, fit.bayes.aif$ktrans, pch=2)
#' abline(0, 1, lwd=2, col=2)
#' legend("right", c("extended/fritz.hansen","orton.exp/user"), pch=1:2)
#' 
#' fit.lm <- dcemri.lm(img, time, mask, aif="fritz.hansen")
#' fit.lm.aif <- dcemri.lm(img, time, mask, model="orton.exp", aif="user",
#'                         user=aifparams)
#' 
#' plot(breast$ktrans, fit.bayes$ktrans, xlim=c(0,1), ylim=c(0,1),
#'      xlab=expression(paste("True ", K^{trans})),
#'      ylab=expression(paste("Estimated ", K^{trans})))
#' points(breast$ktrans, fit.bayes.aif$ktrans, pch=2)
#' points(breast$ktrans, fit.lm$ktrans, pch=3)
#' points(breast$ktrans, fit.lm.aif$ktrans, pch=4)
#' abline(0, 1, lwd=2, col=2)
#' legend("bottomright", c("Bayesian Estimation (fritz-hansen)",
#'                         "Bayesian Estimation (orton.exp)",
#'                         "Levenburg-Marquardt (weinmann/fritz.hansen)",
#'                         "Levenburg-Marquardt (orton.exp/user)"), pch=1:4)
#' 
#' @docType methods
#' @rdname dcemri.bayes-methods
setGeneric("dcemri.bayes", function(conc, ...) standardGeneric("dcemri.bayes"))
#' @export
#' @rdname dcemri.bayes-methods
#' @aliases dcemri.bayes,array-method
#' @useDynLib dcemriS4 dce_bayes_run_single
setMethod("dcemri.bayes", signature(conc="array"),
          function(conc, time, img.mask, model="extended",
                         aif=NULL, user=NULL, nriters=3000, thin=3,
                         burnin=1000, tune=267, ab.ktrans=c(0,1),
                         ab.kep=ab.ktrans, ab.vp=c(1,19),
                         ab.tauepsilon=c(1,1/1000), samples=FALSE,
                         multicore=FALSE, verbose=FALSE, dic=FALSE, ...)
          .dcemriWrapper("dcemri.bayes", conc, time, img.mask, model,
                        aif, user, nriters, thin, burnin, tune, ab.ktrans,
                        ab.kep, ab.vp, ab.tauepsilon, samples,
                        multicore, verbose, dic, ...))

.dcemri.bayes.single <- function(conc, time, nriters=3000, thin=3,
                                 burnin=1000, tune=267, ab.gamma=c(0,1),
                                 ab.theta=c(0,1), ab.vp=c(1,19),
                                 ab.tauepsilon=c(1,1/1000), aif.model=0,
                                 aif.parameter=c(2.4,0.62,3,0.016), vp=1) {

  if (sum(is.na(conc)) > 0) {
    return(NA)
  } else {
    n <- floor((nriters - burnin) / thin)
    if (tune > nriters/2) {
      tune <- floor(nriters/2)
    }
    n0 <- rep(0, n)
    singlerun <- .C("dce_bayes_run_single",
                    as.integer(c(nriters, thin, burnin, tune)),
                    as.double(conc),
                    as.double(ab.gamma),
                    as.double(ab.theta),
                    as.double(ab.vp),
                    as.double(ab.tauepsilon),
                    as.double(c(aif.model, aif.parameter)),
                    as.integer(vp), # is this correct?
                    as.double(time),
                    as.integer(length(time)),
                    as.double(n0),
                    as.double(n0),
                    as.double(n0),
                    as.double(n0),
                    as.double(n0),
                    as.integer(n),
                    PACKAGE="dcemriS4")
    list("ktrans"= singlerun[[11]],
         "kep"= singlerun[[12]],
         "vp"= singlerun[[13]],
         "sigma2"= 1/singlerun[[14]],
         "deviance" = singlerun[[15]])
  }
}

.dcemri.bayes <- function(conc, time, img.mask, model="extended",
                          aif=NULL, user=NULL, nriters=3000, thin=3,
                          burnin=1000, tune=267, ab.ktrans=c(0,1),
                          ab.kep=ab.ktrans, ab.vp=c(1,19),
                          ab.tauepsilon=c(1,1/1000), samples=FALSE,
                          multicore=FALSE, verbose=FALSE, dic=FALSE,
                          ...) {

  ## dcemri.bayes - a function for fitting 1-compartment PK models to
  ## DCE-MRI images using Bayes inference
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
  ##         (ktranserror and keperror), samples if samples=TRUE
  ##

  extractSamples <- function(sample, img.mask, NRI) {
    out <- array(NA, c(NRI, dim(img.mask)))
    out[img.mask] <- sample
    aperm(out, c(2:length(dim(out)), 1)) # not too sure about drop()
  }

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

  if (J > 1 && K > 1) {
    if (sum(dim(img.mask) - dim(conc)[-length(dim(conc))]) != 0) {
      stop("Dimensions of \"conc\" do not agree with \"img.mask\"")
    }
  }
  if (nriters < thin) {
    stop("Please check settings for nriters")
  }
  if (burnin < 0) {
    stop("Please check settings for burnin")
  }
  if (thin < 1) {
    stop("Please check settings for thin")
  }
  if (tune < 50) {
    stop("Please check settings for tune")
  }
  if (burnin < tune) {
    burnin <- tune
    nriters <- nriters + tune
  } else {
    nriters <- nriters + burnin
  }

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
  am <- grep("^[Aa]|^[Mm][^Ee]", names(p))
  aif.parameter <- unlist(p[am])
  if (!is.null(p$D) && p$D != 1) {
    a <- grep("^[Aa]", names(p))
    aif.parameter[a] <- p$D * aif.parameter[a]
  }

  ## translate "model" to "aif.model" and "vp.do"
  switch(model,
         weinmann = {
           aif.model <- 0
           vp.do <- 0
         },
         extended = {
           aif.model <- 0
           vp.do <- 1
         },
         orton.exp = {
           aif.model <- 1
           vp.do <- 1
         },
         kety.orton.exp = {
           aif.model <- 1
           vp.do <- 0
         },
         stop("Model is not supported."))

  ## img.mask <- array(img.mask,c(I,J,K)) # why?
  nvoxels <- sum(img.mask)
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  img.mask <- ifelse(img.mask > 0, TRUE, FALSE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0
  conc.list <- vector("list", nvoxels)
  for (i in 1:nvoxels) {
    conc.list[[i]] <- conc.mat[i,]
  }
  if (verbose) {
    cat("  Estimating the kinetic parameters...", fill=TRUE)
  }
  if (multicore) {
    bayes.list <- parallel::mclapply(conc.list, FUN=.dcemri.bayes.single,
                                     time=time, nriters=nriters, thin=thin,
                                     burnin=burnin, tune=tune, ab.gamma=ab.ktrans,
                                     ab.theta=ab.kep, ab.vp=ab.vp,
                                     ab.tauepsilon=ab.tauepsilon, aif.model=aif.model,
                                     aif.parameter=aif.parameter, vp=vp.do, 
                                     mc.preschedule=TRUE)
  } else {
    bayes.list <- lapply(conc.list, FUN=.dcemri.bayes.single, time=time,
                         nriters=nriters, thin=thin, burnin=burnin,
                         tune=tune, ab.gamma=ab.ktrans, ab.theta=ab.kep,
                         ab.vp=ab.vp, ab.tauepsilon=ab.tauepsilon,
                         aif.model=aif.model, aif.parameter=aif.parameter,
                         vp=vp.do)
  }
  rm(conc.list) ; gc()

  if (verbose) {
    cat("  Extracting results...", fill=TRUE)
  }

  n <- floor((nriters - burnin) / thin)  # number of samples from posterior
  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  sigma2 <- rep(NA, nvoxels)
  if (model %in% c("extended", "orton.exp", "orton.cos")) {
    Vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  }
  if (dic) {
    med.deviance <- rep(NA, nvoxels)
  }
  if (samples) {
    sigma2.samples <- ktrans.samples <- kep.samples <- rep(NA, n*nvoxels)
    if (model %in% c("extended", "orton.exp", "orton.cos")) {
      Vp.samples <- rep(NA, n*nvoxels)
    }
    if (dic) {
      deviance.samples <- rep(NA, n*nvoxels)
    }
  }

  for (k in 1:nvoxels) {
    index <- (k-1)*n + (1:n)
    ktrans$par[k] <- median(bayes.list[[k]]$ktrans)
    kep$par[k] <- median(bayes.list[[k]]$kep)
    ktrans$error[k] <- sd(bayes.list[[k]]$ktrans)
    kep$error[k] <- sd(bayes.list[[k]]$kep)
    if (model %in% c("extended", "orton.exp", "orton.cos")) {
      Vp$par[k] <- median(bayes.list[[k]]$vp)
      Vp$error[k] <- sd(bayes.list[[k]]$vp)
      if (samples) {
        Vp.samples[index] <- bayes.list[[k]]$v
      }
    }
    sigma2[k] <- median(bayes.list[[k]]$sigma2)
    if (dic) {
      med.deviance[k] <- median(bayes.list[[k]]$deviance)
    }
    if (samples) {
      ktrans.samples[index] <- bayes.list[[k]]$ktrans
      kep.samples[index] <- bayes.list[[k]]$kep
      sigma2.samples[index] <- bayes.list[[k]]$sigma2
      if (dic) {
        deviance.samples[index] <- bayes.list[[k]]$deviance
      }
    }
  }
  rm(bayes.list) ; gc()

  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }

  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par
  B[img.mask] <- ktrans$error
  ktrans.out <- list(par=A, error=B)
  rm(A,B,ktrans)
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par
  B[img.mask] <- kep$error
  kep.out <- list(par=A, error=B)
  rm(A,B,kep)
  if (model %in% c("extended", "orton.exp", "orton.cos")) {
    A <- B <- array(NA, c(I,J,K))
    A[img.mask] <- Vp$par
    B[img.mask] <- Vp$error
    Vp.out <- list(par=A, error=B)
    rm(A,B,Vp)
  }

  A <- array(NA, c(I,J,K))
  A[img.mask] <- sigma2
  sigma2.out <- A
  rm(A,sigma2)

  if (dic) {
    A <- array(NA, c(I,J,K))
    A[img.mask] <- med.deviance
    med.deviance <- A
    rm(A)
  }

  if (samples) {
    if (verbose) {
      cat("  Reconstructing samples...", fill=TRUE)
    }
    ktrans.out$samples <- extractSamples(ktrans.samples, img.mask, n)
    kep.out$samples <- extractSamples(kep.samples, img.mask, n)
    if (model %in% c("extended", "orton.exp", "orton.cos")) {
      Vp.out$samples <- extractSamples(Vp.samples, img.mask, n)
    }
    sigma2.samples <- extractSamples(sigma2.samples, img.mask, n)
    if (dic) {
      deviance.samples <- extractSamples(deviance.samples, img.mask, n)
    }
  }

  returnable <- list(ktrans=ktrans.out$par,
                     kep=kep.out$par,
                     ktranserror=ktrans.out$error,
                     keperror=kep.out$error,
                     ve=ktrans.out$par/kep.out$par,
                     time=time)
  if (model %in% c("extended", "orton.exp", "orton.cos")) {
    returnable[["vp"]] <- Vp.out$par
    returnable[["vperror"]] <- Vp.out$error
    if (samples) {
      returnable[["vp.samples"]] <- Vp.out$samples
    }
  }
  ## DIC
  if (dic) {
    if (verbose) {
      cat("  Computing DIC...", fill=TRUE)
    }
    fitted <- array(NA, c(I,J,K,length(time)))
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K) {
          if (img.mask[i,j,k]) {
            par <- NULL
            if (vp.do) {
              par["vp"] <- Vp.out$par[i,j,k]
            }		# Caution: order of assignment is important!!!
            par["ktrans"]=ktrans.out$par[i,j,k]
            par["kep"]=kep.out$par[i,j,k]
            fitted[i,j,k,] <- kineticModel(time, par, model=model, aif=aif)
          }
        }
      }
    }
    conc <- array(conc, c(I,J,K,length(time)))
    fitted <- fitted - conc
    fitted <- fitted * fitted
    fitted <- apply(fitted, 1:3, sum)
    deviance.med <- length(time) * log(sigma2.out) + fitted / sigma2.out
    pD <- med.deviance - deviance.med
    DIC <- med.deviance + pD
    returnable[["DIC"]] <- sum(DIC,na.rm=TRUE)
    returnable[["pD"]] <- sum(pD,na.rm=TRUE)
    returnable[["DIC.map"]] <- DIC
    returnable[["pD.map"]] <- pD
    returnable[["deviance.med"]] <- deviance.med
    returnable[["med.deviance"]] <- med.deviance
    if (samples) {
      returnable[["deviance.samples"]] <- deviance.samples
    }
    rm(DIC)
    rm(pD)
    rm(deviance.med)
    rm(med.deviance)
    rm(deviance.samples)
    gc()
  }
  if (samples) {
    temp <- ktrans.out$samples
    rm(ktrans.out)
    returnable[["ktrans.samples"]] <- temp
    temp <- kep.out$samples
    rm(kep.out)
    returnable[["kep.samples"]] <- temp
    returnable[["sigma2.samples"]] <- sigma2.samples
  }
  ## rm(Vp.out) ; gc()

  return(returnable)
}

