##
## Copyright (c) 2009-2011 Brandon Whitcher and Volker Schmid
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
## $Id: flipangle.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## dam() = double-angle method
#############################################################################
#' Double-Angle Method for B1+ Mapping
#' 
#' For in vivo MRI at high field (\eqn{\geq3}{>=3} T) it is essential to
#' consider the homogeneity of the active B1 field (B1+).  The B1+ field is the
#' transverse, circularly polarized component of B1 that is rotating in the
#' same sense as the magnetization.  When exciting or manipulating large
#' collections of spins, nonuniformity in B1+ results in nonuniform treatment
#' of spins.  This leads to spatially varying image signal and image contrast
#' and to difficulty in image interpretation and image-based quantification.
#' 
#' The proposed method uses an adaptation of the double angle method (DAM).
#' Such methods allow calculation of a flip-angle map, which is an indirect
#' measure of the B1+ field.  Two images are acquired: \eqn{I_1}{I1} with
#' prescribed tip \eqn{\alpha_1}{alpha1} and \eqn{I_2}{I2} with prescribed tip
#' \eqn{\alpha_2=2\alpha_1}{alpha2 = 2*alpha1}.  All other signal-affecting
#' sequence parameters are kept constant. For each voxel, the ratio of
#' magnitude images satisfies
#' \deqn{\frac{I_2(r)}{I_1(r)}=\frac{\sin\alpha_2(r)f_2(T_1,\mbox{TR})}{\sin\alpha_1(r)f_1(T_1,\mbox{TR})}}{}
#' where \eqn{r} represents spatial position and \eqn{alpha_1(r)}{alpha1(r)}
#' and \eqn{\alpha_2(r)}{alpha2(r)} are tip angles that vary with the spatially
#' varying B1+ field.  If the effects of \eqn{T_1}{T1} and \eqn{T_2}{T2}
#' relaxation can be neglected, then the actual tip angles as a function of
#' spatial position satisfy
#' \deqn{\alpha(r)=\mbox{arccos}\left(\left|\frac{I_2(r)}{2I_1(r)}\right|\right)}{}
#' A long repetition time (\eqn{TR\leq{5T_1}}{TR <= 5*T1}) is typically used
#' with the double-angle methods so that there is no \eqn{T_1}{T1} dependence
#' in either \eqn{I_1}{I1} or \eqn{I_2}{I2} (i.e.,
#' \eqn{f_1(T_1,TR)=f_2(T_1,TR)=1.0}{f1(T1,TR) = f2(T1,TR) = 1.0}).  Instead,
#' the proposed method includes a magnetization-reset sequence after each data
#' acquisition with the goal of putting the spin population in the same state
#' regardless of whether the or \eqn{\alpha_2}{alpha2} excitation was used for
#' the preceding acquisition (i.e.,
#' \eqn{f_1(T_1,TR)=f_2(T_1,TR)\ne1.0}{f1(T1,TR) = f2(T1,TR) != 1.0}).
#' 
#' @param low is the (3D) array of signal intensities at the low flip angle.
#' @param high is the (3D) array of signal intensities at the high flip angle
#' (note, 2*low = high).
#' @param low.deg is the low flip angle (in degrees).
#' @return An array, the same dimension as the acquired signal intensities, is
#' returned containing the multiplicative factor associated with the low flip
#' angle acquisition.  That is, if no B1+ inhomogeneity was present then the
#' array would only contain ones.  Numbers other than one indicate the extent
#' of the inhomogeneity as a function of spatial location.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @references 
#' Cunningham, C.H., Pauly, J.M. and Nayak, K.S. (2006) Saturated Double-Angle 
#' Method for Rapid B1+ Mapping, \emph{Magnetic Resonance in Medicine}, 
#' \bold{55}, 1326-1333.
#' @keywords models
#' @export doubleAngleMethod
doubleAngleMethod <- function(low, high, low.deg) {
  alpha <- acos(abs(high /(2*low)))
  (180/pi * alpha) / low.deg # radians to degrees
}

#############################################################################
## R10.lm() = estimate R1 using Levenburg-Marquardt
#############################################################################
#' @export
#' @rdname relaxation-methods
#' @aliases R1.lm
R10.lm <- function(signal, alpha, TR, guess,
                   control=minpack.lm::nls.lm.control()) {
  func <- function(x, signal, alpha, TR) {
    R1 <- x[1]
    m0 <- x[2]
    theta <- pi/180 * alpha # degrees to radians
    signal -
      m0 * sin(theta) * (1 - exp(-TR*R1)) / (1 - cos(theta) * exp(-TR*R1))
  }
  out <- minpack.lm::nls.lm(par=guess, fn=func, control=control,
                            signal=signal, alpha=alpha, TR=TR)
  list(R1=out$par[1], m0=out$par[2], hessian=out$hessian, info=out$info,
       message=out$message)
}

#############################################################################
## E10.lm() = estimate exp(-TR*R1) using Levenburg-Marquardt
#############################################################################
#' @export
#' @rdname relaxation-methods
#' @aliases E10.lm
E10.lm <- function(signal, alpha, guess,
                   control=minpack.lm::nls.lm.control()) {
  func <- function(x, signal, alpha) {
    E1 <- x[1]
    m0 <- x[2]
    theta <- pi/180 * alpha # degrees to radians
    signal - m0 * sin(theta) * (1 - E1) / (1 - cos(theta) * E1)
  }
  out <- minpack.lm::nls.lm(par=guess, fn=func, control=control,
                            signal=signal, alpha=alpha)
  list(E10=out$par[1], m0=out$par[2], hessian=out$hessian, info=out$info,
       message=out$message)
}

#############################################################################
## setGeneric("R1.fast")
#############################################################################
#' Estimate Intrinsic Tissue Relaxivity
#' 
#' Estimation of the intrinsic tissue relaxivity is achieved through nonlinear
#' optimization and the dynamic signal intensities are converted into contrast
#' agent concentration.
#' 
#' The \code{E10.lm} and \code{R10.lm} functions estimate parameters for a
#' vector of observed MR signal intensities, as a function of flip angle, using
#' the following relationship \deqn{S(\alpha) = m_0 \frac{\sin(\alpha) \left(1
#' - \exp{-\textrm{TR}/\textrm{T}_1}\right)}{\left(1 - \cos(\alpha)
#' \exp{-\textrm{TR}/\textrm{T}_1}\right)}.} The only difference between the
#' two functions is exactly what is being estimated in the nonlinear least
#' squares formulation.  They both require the function
#' \code{\link[minpack.lm]{nls.lm}} that uses the Levenberg-Marquardt
#' algorithm.
#' 
#' The \code{CA.fast} function calls on \code{R1.fast} to rearrange the assumed
#' multidimensional (2D or 3D) structure of the multiple flip-angle data into a
#' single matrix to take advantage of internal R functions instead of loops
#' when calling \code{E10.lm}.  Conversion of the dynamic signal intensities to
#' contrast agent concentration is performed via \deqn{[Gd] =
#' \frac{1}{r_1}\left(\frac{1}{\textrm{T}_1} -
#' \frac{1}{\textrm{T}_{10}}\right).}
#' 
#' The \code{CA2.fast} function assumes only two flip angles have been acquired
#' and uses an approximation to the nonlinear relationship between signal
#' intensity and flip angle enable to conversion from signal intensity to
#' contrast agent concentration.
#' 
#' @aliases R10.lm E10.lm R1.fast,array-method CA.fast,array-method
#' CA.fast2,array-method R1.fast CA.fast CA.fast2
#' @param signal is the vector of signal intensities as a function of flip
#' angles.
#' @param ... Additional variables defined by the method.  
#' @param alpha is the vector of flip angles (in degrees).
#' @param TR is the relaxation time (in seconds) used in the acquisition of the
#' MRI data.
#' @param guess is the vector of initial values for the parameters of interest:
#' \eqn{m_0}{M0} and \eqn{R_{10}}{R10}.
#' @param control An optional list of control settings for \code{nls.lm}.  See
#' \code{nls.lm.control} for the names of the settable control values and their
#' effect.
#' @param dynamic a multidimensional array of contrast agent concentrations.
#' The last dimension is assumed to be temporal, while the previous dimenions
#' are assued to be spatial.
#' @param flip.mask,dyn.mask is a (logical) multidimensional array that
#' identifies the voxels to be analyzed.
#' @param dangle is the flip angle used to acquire the dynamic MRI data.
#' @param flip a multidimensional array of contrast agent concentrations.  The
#' last dimension is assumed to be a function of the flip angles, while the
#' previous dimenions are assued to be spatial.
#' @param fangles is the vector of flip angles (in degrees).
#' @param r1 is the spin-lattice relaxivity constant (default = 4.39 for 1.5T).
#' For 3T data it may be necessary to adjust this value.
#' @param multicore is a logical variable (default = \code{FALSE}) that allows
#' parallel processing via \pkg{parallel}.
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @return A list structure is produced with (all or some of the) parameter
#' estimates \item{M0}{Scaling factor between signal intensity and T1.}
#' \item{R10}{Pre-injection tissue relaxation rate (3D array);
#' \eqn{R1_{0}=1/T1_{0}}{R10=1/T10}.} \item{R1t}{Time-varying tissue relaxation
#' rate (4D array); \eqn{R1(t)=1/T1(t)}{R1(t)=1/T1(t)}.} \item{conc}{Contrast
#' agent concentration (4D array).} and information about the convergence of
#' the nonlinear optimization routine.
#' @note The longitudinal relaxivity is set, by default, to
#' \eqn{r_1=4(mM\cdot{s})^{-1}}{r1=4/(mM s)} which is a reasonable value for
#' gadolinium contrast agents at 1.5 Tesla.  Double-check the scanning
#' procedure manual to ensure the correct value is used.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{dcemri.lm}}, \code{\link[minpack.lm]{nls.lm}}
#' @references 
#' Buxton, R.B. (2002) \emph{Introduction to Functional Magnetic
#' Resonance Imaging: Principles & Techniques}, Cambridge University Press:
#' Cambridge, UK.
#' 
#' Li, K.-L., Zhu, X.P., Waterton, J. and Jackson, A. (2000) Improved 3D
#' quantiative mapping of blood volume and endothelial permeability in brain
#' tumors, \emph{Journal of Magnetic Resonance Imaging}, \bold{12}, 347-357.
#' 
#' Li, K.-L., Zhu, X.P., Kamaly-Asl, I.D., Checkley, D.R., Tessier, J.J.L.,
#' Waterton, J.C. and Jackson, A. (2000) Quantification of endothelial
#' permeability, leakage space, and blood volume in brain tumors using combined
#' T1 and T2* contrast-enhanced dynamic MR imaging, \emph{Journal of Magnetic
#' Resonance Imaging}, \bold{11}, 575-585.
#' 
#' Parker, G.J.M. and Padhani, A.R. (2003) \eqn{T_1}{T1}-w DCE-MRI:
#' \eqn{T_1}{T1}-weighted Dynamic Contrast-enhanced MRI, in \emph{Quantiative
#' MRI of the Brain} (P. Tofts ed.), Wiley: Chichester, UK, pp. 341-364.
#' @keywords misc
#' @examples
#' 
#' ## Parameters for simulated data
#' S0 <- 100
#' TR <- 5 / 1000            # seconds
#' T1 <- 1.5                 # seconds
#' alpha <- seq(2, 24, by=2) # degrees
#' 
#' ## Signal intensities for spoiled gradient echo image
#' gre <- function(S0, TR, T1, alpha) {
#'   theta <- alpha * pi/180 # radians
#'   S0 * (1 - exp(-TR/T1)) * sin(theta) / (1 - cos(theta) * exp(-TR/T1))
#' }
#' set.seed(1234)
#' signal <- array(gre(S0, TR, T1, alpha) + rnorm(length(alpha), sd=.15),
#'                 c(rep(1,3), length(alpha)))
#' out <- R1.fast(signal, array(TRUE, rep(1,3)), alpha, TR)
#' unlist(out)
#' plot(alpha, signal, xlab="Flip angle", ylab="Signal intensity")
#' lines(alpha, gre(S0, TR, T1, alpha), lwd=2, col=1)
#' lines(alpha, gre(c(out$M0), TR, 1/c(out$R10), alpha), lwd=2, col=2)
#' legend("topright", c("True","Estimated"), lwd=2, col=1:2)
#' 
#' @export
#' @docType methods
#' @rdname relaxation-methods
setGeneric("R1.fast", function(flip, ...) standardGeneric("R1.fast"))
#' @export
#' @rdname relaxation-methods
#' @aliases R1.fast,array-method
setMethod("R1.fast", signature(flip="array"),
          function(flip, flip.mask, fangles, TR,
                   control=minpack.lm::nls.lm.control(),
                   multicore=FALSE, verbose=FALSE)
	    .dcemriWrapper("R1.fast", flip, flip.mask, fangles, TR, control,
                           multicore, verbose))

#############################################################################
## R1.fast()
#############################################################################

.R1.fast <- function(flip, flip.mask, fangles, TR,
                     control=minpack.lm::nls.lm.control(),
                     multicore=FALSE, verbose=FALSE) {

  if (length(dim(flip)) != 4) { # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  }
  if (!is.logical(flip.mask)) { # Check flip.mask is logical
    stop("Mask must be logical.")
  }
  if (length(fangles) > 1 && length(TR) == 1) {
    method <- "E10"
  } else {
    if (length(fangles) == 1 && length(TR) > 1) {
      method <- "R10"
    } else {
      stop("Only vector of flip angles and single TR or single flip angle and vector TR is allowed.")
    }
  }

  X <- nrow(flip)
  Y <- ncol(flip)
  Z <- oro.nifti::nsli(flip)
  nvoxels <- sum(flip.mask)

  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  flip.mat <- matrix(flip[flip.mask], nrow=nvoxels)
  if (is.array(fangles)) {
    fangles.mat <- matrix(fangles[flip.mask], nrow=nvoxels)
  } else {
    fangles.mat <- matrix(fangles, nrow=nvoxels, ncol=length(fangles),
                          byrow=TRUE)
  }
  flip.list <- vector("list", nvoxels)
  for (k in 1:nvoxels) {
    flip.list[[k]] <- list(signal=flip.mat[k,], angles=fangles.mat[k,])
  }
  if (verbose) {
    cat("  Calculating R10 and M0...", fill=TRUE)
  }
  if (multicore) {
    if (method == "E10") {
      fit.list <- parallel::mclapply(flip.list, function(x) {
        E10.lm(x$signal, x$angles, guess=c(1, mean(x$signal)), control)
      })
    } else {
      fit.list <- parallel::mclapply(flip.list, function(x) {
        R10.lm(x$signal, x$angles, TR, guess=c(1, mean(x$signal)), control)
      })
    }
  } else {
    if (method == "E10") {
      fit.list <- lapply(flip.list, function(x) {
        E10.lm(x$signal, x$angles, guess=c(1, mean(x$signal)), control)
      })
    } else {
      fit.list <- lapply(flip.list, function(x) {
        R10.lm(x$signal, x$angles, TR, guess=c(1, mean(x$signal)), control)
      })
    }
  }
  rm(flip.list) ; gc()
  R10 <- M0 <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  for (k in 1:nvoxels) {
    if (fit.list[[k]]$info > 0 && fit.list[[k]]$info < 5) {
      if (method == "E10") {
        R10$par[k] <- log(fit.list[[k]]$E10) / -TR
      } else {
        R10$par[k] <- fit.list[[k]]$R10
      }
      M0$par[k] <- fit.list[[k]]$m0
      R10$error[k] <- sqrt(fit.list[[k]]$hessian[1,1])
      M0$error[k] <- sqrt(fit.list[[k]]$hessian[2,2])
    } else {
      R10$par[k] <- M0$par[k] <- R10$error[k] <- M0$error[k] <- NA
    }
  }

  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  R10.array <- M0.array <- R10error <- M0error <- array(NA, c(X,Y,Z))
  R10.array[flip.mask] <- R10$par
  M0.array[flip.mask] <- M0$par
  R10error[flip.mask] <- R10$error
  M0error[flip.mask] <- M0$error

  list(M0 = M0.array, R10 = R10.array, M0.error = NULL, R10.error = NULL)
}

#############################################################################
## setGeneric("CA.fast")
#############################################################################
#' @export
#' @rdname relaxation-methods
setGeneric("CA.fast", function(dynamic, ...) standardGeneric("CA.fast"))
#' @export
#' @rdname relaxation-methods
#' @aliases CA.fast,array-method
setMethod("CA.fast", signature(dynamic="array"),
	  function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4,
                   control=minpack.lm::nls.lm.control(maxiter=200),
                   multicore=FALSE, verbose=FALSE)
	    .dcemriWrapper("CA.fast", dynamic, dyn.mask, dangle, flip,
                           fangles, TR, r1, control, multicore, verbose))

#############################################################################
## CA.fast() = estimate contrast-agent concentration and other stuff
#############################################################################

.CA.fast <- function(dynamic, dyn.mask, dangle, flip, fangles, TR,
                     r1=4, control=minpack.lm::nls.lm.control(maxiter=200),
                     multicore=FALSE, verbose=FALSE) {

  if (length(dim(flip)) != 4) { # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  }
  if (!is.logical(dyn.mask)) { # Check dyn.mask is logical
    stop("Mask must be logical.")
  }
  if (! identical(length(fangles), oro.nifti::ntim(flip)) &&
      ! isTRUE(all.equal(dim(flip), dim(fangles)))) {
    ## Check that #(flip angles) are equal
    stop("Number of flip angles must agree with dimension of flip-angle data.")
  }

  R1est <- R1.fast(flip, dyn.mask, fangles, TR, control, multicore, verbose)

  if (verbose) {
    cat("  Calculating concentration...", fill=TRUE)
  }
  theta <- dangle * pi/180
  A <- sweep(sweep(dynamic, 1:3, dynamic[,,,1], "-"),
             1:3, R1est$M0, "/") / sin(theta)
  B <- (1 - exp(-TR * R1est$R10)) / (1 - cos(theta) * exp(-TR * R1est$R10))
  AB <- sweep(A, 1:3, B, "+")
  rm(A,B)
  R1t <- -(1/TR) * log((1 - AB) / (1 - cos(theta) * AB))
  rm(AB)
  conc <- sweep(R1t, 1:3, R1est$R10, "-") / r1

  list(M0 = R1est$M0, R10 = R1est$R10, R1t = R1t, conc = conc)
}

#############################################################################
## setGeneric("CA.fast2")
#############################################################################
#' @export
#' @rdname relaxation-methods
setGeneric("CA.fast2", function(dynamic, ...) standardGeneric("CA.fast2"))
#' @export
#' @rdname relaxation-methods
#' @aliases CA.fast2,array-method
setMethod("CA.fast2", signature(dynamic="array"),
	  function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4,
                   verbose=FALSE)
          .dcemriWrapper("CA.fast2", dynamic, dyn.mask, dangle, flip,
                         fangles, TR, r1, verbose))

#############################################################################
## CA.fast2()
#############################################################################

.CA.fast2 <- function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4,
                     verbose=FALSE) {

  if (length(dim(flip)) != 4) {  # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  }
  if (! identical(length(fangles), oro.nifti::ntim(flip)) &&
      ! isTRUE(all.equal(dim(flip), dim(fangles)))) {
    ## Check that #(flip angles) are equal
    stop("Number of flip angles must agree with dimension of flip-angle data.")
  }
  ##if (ntim(flip) != 2 || length(fangles) != 2) {
  ##  stop("Only two flip angles are allowed.")
  ##}
  if (!is.logical(dyn.mask)) { # Check dyn.mask is logical
    stop("Mask must be logical.")
  }
  nangles <- length(fangles)
  nvoxels <- sum(dyn.mask)
  M <- nrow(flip)
  N <- ncol(flip)
  Z <- oro.nifti::nsli(dynamic)
  W <- oro.nifti::ntim(dynamic)
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  if (is.array(fangles)) {
    fangles.mat <- matrix(fangles[dyn.mask], nrow=nvoxels)
  } else {
    fangles.mat <- matrix(fangles, nrow=nvoxels, ncol=length(fangles),
                          byrow=TRUE)
  }
  dyn.mat <- matrix(dynamic[dyn.mask], nvoxels)
  flip.mat <- matrix(flip[dyn.mask], nvoxels)
  R10 <- M0 <- numeric(nvoxels)
  if (verbose) {
    cat("  Calculating R10 and M0...", fill=TRUE)
  }
  x <- flip.mat / tan(pi * fangles.mat / 180)
  x <- ifelse(is.finite(x), x, 0)
  y <- flip.mat / sin(pi * fangles.mat / 180)
  y <- ifelse(is.finite(y), y, 0)
  for (k in 1:nvoxels) {
    #x <- c(flip.mat[k,1] / tan(pi * fangles.mat[k,1] / 180),
    #       flip.mat[k,2] / tan(pi * fangles.mat[k,2] / 180))
    #x <- ifelse(is.finite(x), x, 0)
    #y <- c(flip.mat[k,1] / sin(pi * fangles.mat[k,1] / 180),
    #       flip.mat[k,2] / sin(pi * fangles.mat[k,2] / 180))
    #y <- ifelse(is.finite(y), y, 0)
    fit <- lsfit(x[k, ], y[k, ])$coefficients
    R10[k] <- log(fit[2]) / -TR
    M0[k] <- fit[1] / (1 - fit[2])
  }
  if (verbose) {
    cat("  Calculating concentration...", fill=TRUE)
  }
  theta <- dangle * pi/180
  CD <- conc <- matrix(NA, nvoxels, W)
  B <- (1 - exp(-TR * R10)) / (1 - cos(theta) * exp(-TR * R10))
  A <- (dyn.mat - dyn.mat[,1]) / M0 / sin(theta)
  R1t <- -(1/TR) * log((1 - (A+B)) / (1 - cos(theta) * (A+B)))
  conc <- (R1t - R10) / r1
  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  R10.array <- M0.array <- array(NA, c(M,N,Z))
  R10.array[dyn.mask] <- R10
  M0.array[dyn.mask] <- M0
  conc.array <- R1t.array <- array(NA, c(M,N,Z,W))
  mask4D <- array(dyn.mask, c(M,N,Z,W))
  conc.array[mask4D] <- unlist(conc)
  R1t.array[mask4D] <- unlist(R1t)
  list(M0 = M0.array, R10 = R10.array, R1t = R1t.array, conc = conc.array)
}
