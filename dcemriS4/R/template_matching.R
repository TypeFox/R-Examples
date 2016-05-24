##
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
## $Id:$
##

#' Convolution of 3D Arrays using the Fourier Transform
#' 
#' Convolve a three-dimensinal array with another three-dimensional arry using 
#' the Fast Fourier Transform (FFT).
#' 
#' 
#' The arrays \eqn{A} and \eqn{B} are transformed into the Fourier domain and 
#' multiplied together (equivalent to a convolution in the image domain across 
#' all spatial locations simultaneously).
#' 
#' @param A is a three-dimensional array (\dQuote{the template}).
#' @param B is a three-dimensional array (\dQuote{the target}).
#' @param C is a vector of length three (the center of \dQuote{the template}).
#' @param FFTA is the three-dimensional Fourier transform of \code{A}, this may 
#'   save time when looping over multiple \dQuote{targets}.
#' @return A three-dimensional array, the same dimension as the input arrays, 
#'   that is the convolution of the \dQuote{target} to the \dQuote{template} at 
#'   all spatial locations.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{fft}}, \code{\link{fastTemplateMatching}}, \code{\link{shift3D}}
#' @references 
#' Briggs, W.L. and Henson, V.E. (1995) \emph{The DFT: An Owner's
#' Manual for the Discrete Fourier Transform}, SIAM: Philadelphia.
#' @examples
#' 
#' cube <- array(0, c(20,20,1))
#' cube[9:12,9:12,1] <- 1
#' tkernel <- array(0, c(20,20,1))
#' tkernel[,,1] <- c(.5, 1, .5, rep(0,20-3)) %o% c(.5, 1, .5, rep(0,20-3))
#' tcenter <- findCenter(ifelse(tkernel > 0, TRUE, FALSE))
#' out <- convFFT(tkernel, cube, tcenter)
#' out[8:13,8:13,1]  ## text output
#' 
#' par(mfrow=c(2,2))  ## graphic output
#' image(drop(tkernel), col=oro.nifti::tim.colors(), main="Template")
#' image(drop(cube), col=oro.nifti::tim.colors(), main="Target")
#' image(drop(out), col=oro.nifti::tim.colors(), main="Output")
#' 
#' @export convFFT
convFFT <- function(A, B, C, FFTA=NULL) {
  if (length(dim(A)) == 3) {
    if (length(dim(A)) == length(dim(B)) && length(dim(A)) == length(C)) {
      X <- nrow(A)
      Y <- ncol(A)
      Z <- oro.nifti::nsli(A)
      if (is.null(FFTA)) {
        out <- Re(fft(fft(A) * Conj(fft(B)), inverse=TRUE))[X:1,Y:1,Z:1,drop=FALSE] / (X*Y*Z)
      } else {
        out <- Re(fft(FFTA * Conj(fft(B)), inverse=TRUE))[X:1,Y:1,Z:1,drop=FALSE] / (X*Y*Z)
      }
      if (X > 1) {
        out <- out[c((X-C[1]+1):X, 1:(X-C[1])),,,drop=FALSE]
      }
      if (Y > 1) {
        out <- out[,c((Y-C[2]+1):Y, 1:(Y-C[2])),,drop=FALSE]
      }
      if (Z > 1) {
        out <- out[,,c((Z-C[3]+1):Z, 1:(Z-C[3])),drop=FALSE]
      }
      return(out)
    } else {
      stop("Objects are not all the same dimension!")
    }
  } else {
    stop("Only three-dimensional objects are allowed!")
  }
}

#' Find the Center of a Binary Mask
#' 
#' The center of a binary mask is determined.
#' 
#' This method most likely only works with convex three-dimensional shapes
#' (e.g., a hyper-rectangle).  Further testing is required to know the limits
#' of the current implementation.
#' 
#' @param M is a binary mask (multidimensional array of logical values).
#' @return A vector of values the same length as the input array.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{fastTemplateMatching}}
#' @keywords misc
#' @examples
#' 
#' M <- array(FALSE, rep(10,3))
#' M[6:10,6:10,6:10] <- TRUE
#' Mc <- findCenter(M)
#' print(Mc)
#' 
#' @export findCenter
findCenter <- function(M) {
  if (!is.logical(M)) {
    stop("Object must be logical!")
  }
  if (length(dim(M)) != 3) {
    stop("Object must be three-dimensional!")
  }
  X <- nrow(M)
  Y <- ncol(M)
  Z <- oro.nifti::nsli(M)
  xx <- array(1:X, dim(M))
  yy <- array(rep(1:Y, each=X), dim(M))
  zz <- array(rep(1:Z, each=X*Y), dim(M))
  center <- c(mean(xx[M]), mean(yy[M]), mean(zz[M]))
  trunc(center)
}

### FIXME Should this be genericised? 

#' Shift a 3D Array in One Dimension
#' 
#' One axis of the three-dimensional array is translated by an integer amount.
#' This is useful when applying convolution operators in the Fourier domain.
#' 
#' @param A is a three-dimensional array.
#' @param s is the integer number of translation steps.
#' @param type is a character string using anatomical coordinates assuming a
#' transverse acquisition scheme (\dQuote{LR} = left-right = x-axis,
#' \dQuote{AP} = anterior-posterior = y-axis, \dQuote{SI} = superior-inferior =
#' z-axis).
#' @param fill is the quantity used to fill gaps induced by the translations
#' (circular boundary conditions are NOT used).
#' @return A three-dimensional array is returned, the same dimension as the
#' original array, with one dimension translated.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{convFFT}}
#' @examples
#' 
#' cube <- array(0, rep(20,3))
#' cube[9:12,9:12,9:12] <- 1
#' cube.shift <- shift3D(cube, 5, type="AP")
#' par(mfrow=c(1,2), mar=rep(0.5,4))
#' image(cube[,,10], xlab="", ylab="", axes=FALSE)
#' image(cube.shift[,,10], xlab="", ylab="", axes=FALSE)
#' 
#' @export shift3D
shift3D <- function(A, s, type, fill=0) {
  if (length(dim(A)) != 3) {
    stop("Object must be three-dimensional!")
  }
  X <- nrow(A)
  Y <- ncol(A)
  Z <- oro.nifti::nsli(A)
  if (s != 0) {
    if (type == "LR") {
      ## left-right
      if (s > 0) {
        A <- A[c((X-s+1):X,1:(X-s)),,]
        A[1:s,,] <- fill
      } else {
        A <- A[c((abs(s)+1):X,1:abs(s)),,]
        A[(X-abs(s)+1):X,,] <- fill
      }
    } else {
      if (type == "AP") {
        ## anterior-posterior
        if (s > 0) {
          A <- A[,c((Y-s+1):Y,1:(Y-s)),]
          A[,1:s,] <- fill
        } else {
          A <- A[,c((abs(s)+1):Y,1:abs(s)),]
          A[,(Y-abs(s)+1):Y,] <- fill
        }
      } else {
        if (type == "SI") {
          ## superior-inferior
          if (s > 0) {
            A <- A[,,c((Z-s+1):Z,1:(Z-s))]
            A[,,1:s] <- 0
          } else {
            A <- A[,,c((abs(s)+1):Z,1:abs(s))]
            A[,,(Z-abs(s)+1):Z] <- fill
          }
        } else {
          stop("Type of translation not recognized!")
        }
      }
    }
  }
  return(A)
}

#############################################################################
## setGeneric("fastTemplateMatching")
#############################################################################
#' Fast Template Matching via Cross-Correlation
#' 
#' Motion correction and/or co-registration of three-dimensional arrays
#' (medical imaging data) are performed by applying a user-defined mask of
#' voxels.  Normalized cross-correlations (in 3D) are computed using the FFT.
#' 
#' An extremely basic method of motion correction/co-registration is
#' implemented by estimating \dQuote{local} cross-correlations based on a
#' binary mask that is a subset of the original three-dimensional volume.  All
#' convolutions are preformed via the FFT (\code{\link{fft}}) and repetitive
#' calculations are minimized where possible.
#' 
#' Only whole-voxel translations are considered.  This does not begin to
#' capture the true effects of motion in soft tissue, but we assume that the
#' object of interest (e.g., tumor) is a fairly rigid structure.  Potential
#' extensions include rigid-body, affine and nonlinear registration techniques
#' along with interploation schemes in order to capture intra-voxel
#' manipulations of the data.
#' 
#' @aliases fastTemplateMatching fastTemplateMatching,array-method
#' @param input is a four-dimensional array of signal intensities.
#' @param ... Additional variables passed to the \code{plot} function.
#' @return A list of objects are returned: \item{out}{Motion-corrected version
#' of the four-dimensional array.} \item{offset}{Translations (in 3D) for each
#' volume in the 4D array.} \item{t.center}{Estimated center of the binary
#' mask.}
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{convFFT}}, \code{\link{findCenter}},
#' \code{\link{shift3D}}
#' @references 
#' Lewis, J.P. (2003) Fast normalized cross-correlation.\cr
#' \url{www.idiom.com/~zilla/}
#' @rdname fastTemplateMatching
#' @docType methods
#' @export
setGeneric("fastTemplateMatching", 
           function(input, ...) standardGeneric("fastTemplateMatching"))
#' @rdname fastTemplateMatching
#' @export
setMethod("fastTemplateMatching", signature(input="array"),
          function(input, ...) .dcemriWrapper("fastTemplateMatching", input, ...))

.ftm <- function(input, mask, reference, plot=TRUE, ...) {
  ## Fast template matching via cross-correlation
  W <- oro.nifti::ntim(input)
  if (W < 1) {
    stop("4D object is assumed!")
  }
  template <- ifelse(mask > 0, reference, 0)
  templateFFT <- fft(template)
  maskFFT <- fft(mask)
  tc <- findCenter(ifelse(template > 0, TRUE, FALSE))
  templateSS <- convFFT(mask, template^2, tc, FFTA=maskFFT)
  ##numerator <- localSS <- array(0, dim(input))
  localCOR <- array(0, dim(input))
  for (w in 1:W) {
    target <- input[,,,w]
    ##numerator[,,,w] <- convFFT(template, target, tc, FFTA=templateFFT)
    ##localSS[,,,w] <- convFFT(mask, target^2, tc, FFTA=maskFFT)
    numerator <- convFFT(template, target, tc, FFTA=templateFFT)
    localSS <- convFFT(mask, target*target, tc, FFTA=maskFFT)
    localCOR[,,,w] <- numerator / sqrt(templateSS[tc[1],tc[2],tc[3]]) / sqrt(localSS)
  }
  ## localCOR <- numerator / sqrt(templateSS[tc[1],tc[2],tc[3]]) / sqrt(localSS)
  localCOR[!is.finite(localCOR)] <- NA
  wmax.localCOR <- apply(localCOR, 4,
                         function(x) {
                           which(x == max(x,na.rm=TRUE), arr.ind=TRUE)
                         })
  ## which(localCOR[,,,w] == max(localCOR[,,,w], na.rm=TRUE), arr.ind=TRUE)
  offset <- tc - wmax.localCOR
  output <- input@.Data
  for (w in 1:W) {
    output[,,,w] <- shift3D(output[,,,w], offset[1,w], type="LR")
    output[,,,w] <- shift3D(output[,,,w], offset[2,w], type="AP")
    output[,,,w] <- shift3D(output[,,,w], offset[3,w], type="SI")
  }
  storage.mode(output) <- storage.mode(input)
  as(output, "nifti") <- input
  if (plot) {
    matplot(1:ncol(offset), t(offset), type="l", lwd=2, lty=1,
            ylim=range(c(range(offset),-10,10)), xlab="", ylab="voxels",
            main="Motion Correction via Template Matching")
    legend("topleft", c("X","Y","Z"), col=1:3, lty=1, lwd=2, bty="n")
  }
  list(out=output, offset=offset, t.center=tc, localCOR=localCOR, plot=NA)
}
