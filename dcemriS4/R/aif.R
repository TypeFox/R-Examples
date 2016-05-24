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
##
## Time-stamp: <2009-07-14 09:47:30 (bjw34032)>
## $Id: aif.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#' Arterial Input Functions
#' 
#' Parametric models for arterial input functions (AIFs) that are compatible
#' with single compartment models for dynamic contrast-enhanced MRI (DCE-MRI).
#' 
#' \code{aif.orton.exp} displays the exponential AIF from Orton \emph{et al.}
#' (2008) for a known set of AIF parameter values.  \code{model.orton.exp}
#' displays the exponential AIF from Orton \emph{et al.} (2008) for a known set
#' of AIF and compartmental model parameter values.  \code{orton.exp.lm}
#' estimates the AIF parameters, using nonlinear optimization, using a vector
#' of observed contrast agent concentrations.
#' 
#' @name aif-models
#' @aliases aif.orton.exp model.orton.exp orton.exp.lm
#' @param tt is a vector of acquisition times (in minutes) relative to
#' injection of the contrast agent.  Negative values should be used prior to
#' the injection.
#' @param AB,muB,AG,muG are parameters of the double exponential function that
#' describe the AIF.
#' @param aparams is the vector of parameters (\eqn{A_B}, \eqn{\mu_B},
#' \eqn{A_G}, \eqn{\mu_G}) associated with the AIF.
#' @param kparams is the vector of parameters (\eqn{v_p}, \eqn{K^{trans}},
#' \eqn{k_{ep}}) associated with the \dQuote{extended Kety model} for contrast
#' agent concentration.
#' @param aif is the vector of observed contrast agent concentrations (data)
#' used to estimate the parametric model.
#' @param guess Initial parameter values for the nonlinear optimization.
#' @param nprint is an integer, that enables controlled printing of iterates if
#' it is positive.  In this case, estimates of \code{par} are printed at the
#' beginning of the first iteration and every \code{nprint} iterations
#' thereafter and immediately prior to return.  If \code{nprint} is not
#' positive, no tracing information on the progress of the optimization is
#' produced.
#' @return \code{aif.orton.exp} and \code{model.orton.exp} return the AIF
#' associated with the pre-specified parameter values.
#' 
#' \code{orton.exp.lm} returns a list structure with \item{AB}{The amplitude of
#' the first exponential function.} \item{muB}{The decay rate of the first
#' exponential function.} \item{AG}{The amplitude of the second exponential
#' function.} \item{muG}{The decay rate of the second exponential function.}
#' \item{info}{The success (or failure) code from the Levenburg-Marquardt
#' algorithm \code{nls.lm}.} \item{message}{The text message associated with
#' the \code{info} paramters.}
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{dcemri.lm}}, \code{\link{extractAIF}},
#' \code{\link[minpack.lm]{nls.lm}}
#' @references Orton, M.R., Collins, D.J., Walker-Samuel, S., d'Arcy, J.A.,
#' Hawkes, D.J., Atkinson, D. and Leach, M.O. (2007) Bayesian estimation of
#' pharmacokinetic parameters for DCE-MRI with a robust treatment of
#' enhancement onset time, \emph{Physics in Medicine and Biology} \bold{52},
#' 2393-2408.
#' 
#' Orton, M.R., d'Arcy, J.A., Walker-Samuel, S., Hawkes, D.J., Atkinson, D.,
#' Collins, D.J. and Leach, M.O. (2008) Computationally efficient vascular
#' input function models for quantitative kinetic modelling using DCE-MRI,
#' \emph{Physics in Medicine and Biology} \bold{53}, 1225-1239.
#' @keywords models
#' @examples
#' 
#' data("buckley")
#' ## Generate AIF params using the orton.exp function from Buckley's AIF
#' xi <- seq(5, 300, by=5)
#' time <- buckley$time.min[xi]
#' aif <- buckley$input[xi]
#' aifparams <- orton.exp.lm(time, aif)
#' aifparams$D <- 1 
#' unlist(aifparams[1:4])
#' 
#' aoe <- aif.orton.exp(time, aifparams$AB, aifparams$muB, aifparams$AG,
#'                      aifparams$muG)
#' with(buckley, plot(time.min, input, type="l", lwd=2))
#' lines(time, aoe, lwd=2, col=2)
#' legend("right", c("Buckley's AIF", "Our approximation"), lty=1,
#'        lwd=2, col=1:2)
#' cbind(time, aif, aoe)[1:10,]
#' 
#' @rdname aif-models
#' @export 
aif.orton.exp <- function(tt, AB, muB, AG, muG) {
  out <- AB * tt * exp(-muB * tt) + AG * (exp(-muG * tt) - exp(-muB * tt))
  out[tt < 0] <- 0
  return(out)
}
#' @rdname aif-models
#' @export orton.exp.lm
orton.exp.lm <- function(tt, aif,
                         guess=c(log(100), log(10), log(1), log(0.1)),
                         nprint=0) {
  func <- function(x, aparams, aif) {
    AB <- aparams[1]
    muB <- aparams[2]
    AG <- aparams[3]
    muG <- aparams[4]
    return(aif - aif.orton.exp(x, AB, muB, AG, muG))
  }
  out <- minpack.lm::nls.lm(par=guess, fn=func, control=list(nprint=nprint),
                            x=tt, aif=aif)
  list(AB=out$par[1], muB=out$par[2], AG=out$par[3], muG=out$par[4], 
       info=out$info, message=out$message)
}
#' @rdname aif-models
#' @export model.orton.exp
model.orton.exp <- function(tt, aparams, kparams) {
  ## Extended model using the exponential AIF from Matthew Orton (ICR)
  Cp <- function(tt, ...) {
    AB * tt * exp(-muB * tt) + AG * (exp(-muG * tt) - exp(-muB * tt))
  }
  aparams <- as.numeric(aparams)
  kparams <- as.numeric(kparams)
  AB <- aparams[1]
  muB <- aparams[2]
  AG <- aparams[3]
  muG <- aparams[4]
  vp <- kparams[1]
  ktrans <- kparams[2]
  kep <- kparams[3]
  
  T1 <- AB * kep / (kep - muB)
  T2 <- tt * exp(-muB * tt) -
    (exp(-muB * tt) - exp(-kep * tt)) / (kep - muB)
  T3 <- AG * kep
  T4 <- (exp(-muG * tt) - exp(-kep * tt)) / (kep - muG) -
    (exp(-muB * tt) - exp(-kep * tt)) / (kep - muB)
  
  out <- vp * Cp(tt) + ktrans * (T1 * T2 + T3 * T4)
  out[tt <= 0] <- 0
  return(out)
}
#' Seed Growing for a 4D Array
#' 
#' Seed growing algorithm to find voxels in a three-dimensional array according 
#' to their correlation to a seed voxel.  The correlation is measured according 
#' to the fourth dimension of the array.
#' 
#' Correlation coefficients are computed for every voxel in the input array.  A 
#' recursive algorithm is then used to grow the region of interest 
#' (\acronym{ROI}) from the seed voxel in three dimensions.  All adjacent 
#' voxels, where the correlation exceeds the threshold, are included.
#' 
#' @aliases extractAIF
#' @param img is the four-dimensional array of medical imaging data.
#' @param x,y,z are the coordinates of the seed voxel.
#' @param thresh is the minimum correlation for inclusion in the region.
#' @return 
#' \item{coord}{is a matrix of the three-dimesional coordinates
#' \eqn{(x,y,z)} for all voxels found by the algorithm.}
#' \item{conc}{is a matrix whose rows correspond to the voxels found by the 
#' algorithm and whose columns are the fourth dimension from the input array 
#' (e.g., contrast agent concentration time curve).}
#' \item{mask}{is an array of boolean values, where only voxels included by the
#' algorithm are given a value greater than zero.} 
#' \item{cor}{is an array that mimics the \code{mask}, but contains the 
#' estimated correlation coefficients for all voxels included by the algorithm.}
#' @author Volker Schmid \email{volker.schmid@@users.sourceforge.net}
#' @keywords misc
#' @rdname extractAIF
#' @export
extractAIF <- function(img, x, y, z, thresh=0.9) {
  c.start <- function(ctc) {
    if (sum(is.na(ctc)) > 0) {
      return(0)
    } else {
      if (sd(ctc) == 0) {
        return(0)
      } else {
        out <- cor(ctc, start, use="pairwise.complete.obs")
        return(out)
      }
    }
  }
  
  check <- function(xx, yy, zz, aif.mask, thresh) {
    if (xx != 1 && aif.mask[xx-1,yy,zz] == 0) {
      if (c.test[xx-1,yy,zz] > thresh) {
        aif.mask[xx-1,yy,zz] <- 1
        aif.mask <- check(xx-1, yy, zz, aif.mask, thresh)
      }
    }
    if (xx != X && aif.mask[xx+1,yy,zz] == 0) {
      if (c.test[xx+1,yy,zz] > thresh) {
        aif.mask[xx+1,yy,zz] <- 1
        aif.mask <- check(xx+1, yy, zz, aif.mask, thresh)
      }
    }
    if (yy != 1 && aif.mask[xx,yy-1,zz] == 0) {
      if (c.test[xx,yy-1,zz] > thresh) {
        aif.mask[xx,yy-1,zz] <- 1
        aif.mask <- check(xx, yy-1, zz, aif.mask, thresh)
      }
    }
    if (yy != Y && aif.mask[xx,yy+1,zz] == 0) {
      if (c.test[xx,yy+1,zz] > thresh) {
        aif.mask[xx,yy+1,zz] <- 1
        aif.mask <- check(xx, yy+1, zz, aif.mask, thresh)
      }
    }
    if (zz != 1 && aif.mask[xx,yy,zz-1] == 0) {
      if (c.test[xx,yy,zz-1] > thresh) {
        aif.mask[xx,yy,zz-1] <- 1
        aif.mask <- check(xx, yy, zz-1, aif.mask, thresh)
      }
    }
    if (zz != Z && aif.mask[xx,yy,zz+1] == 0) {
      if (c.test[xx,yy,zz+1] > thresh) {
        aif.mask[xx,yy,zz+1] <- 1
        aif.mask <- check(xx, yy, zz+1, aif.mask, thresh)
      }
    }
    return(aif.mask)
  }

  X <- dim(img)[1]
  Y <- dim(img)[2]
  Z <- dim(img)[3]
  W <- dim(img)[4]
  
  start <- img[x,y,z,]
  c.test <- apply(img, 1:3, c.start)
  
  aif.mask <- array(0, c(X,Y,Z))
  aif.mask[x,y,z] <- 1
  
  aif.mask <- check(x, y, z, aif.mask, thresh)
  n <- sum(aif.mask, na.rm=TRUE)
  test <- array(NA, c(n,W))
  coord <- array(NA, c(n,3))
  l <- 0
  for (i in 1:X) {
    for (j in 1:Y) {
      for (k in 1:Z) {
        if (!is.na(c.test[i,j,k]) && aif.mask[i,j,k] == 1) {
          l <- l + 1
          test[l,] <- img[i,j,k,]
          coord[l,] <- c(i,j,k)
        }
      }
    }
  }
  if (l != n) {
    return(FALSE)
  }
  list("coord"=coord, "conc"=test, "mask"=aif.mask, "cor"=c.test)
}

