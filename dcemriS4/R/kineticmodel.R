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
## $Id: kineticmodel.R 332 2010-01-29 16:54:07Z bjw34032 $

#' Pharmacokinetic Models
#' 
#' Kinetic curves from single compartment models are computed from kinetic
#' parameters.
#' 
#' Compartmental models are the solution to the modified general rate equation
#' (Kety 1951).  The specific parametric models considered here include the
#' basic Kety model
#' \deqn{C_t(t)=K^{trans}\left[C_p(t)\otimes\exp(-k_{ep}t)\right],} where
#' \eqn{\otimes}{o} is the convolution operator, and the so-called extended
#' Kety model
#' \deqn{C_t(t)=v_pC_p(t)+K^{trans}\left[C_p(t)\otimes\exp(-k_{ep}t)\right].}
#' The arterial input function must be literature-based (with fixed
#' parameters).
#' 
#' @param time is a vector of acquisition times (in minutes).
#' @param par is a list of kinetic parameters; e.g.,
#' \code{list("ktrans"=0.5,"kep"=1)}.
#' @param model is a character string that identifies the type of compartmental
#' model to be used.  Acceptable models include: \dQuote{weinmann} Tofts &
#' Kermode AIF convolved with single compartment model \dQuote{extended}
#' (default) Weinmann model extended with additional vascular compartment, ...
#' @param aif is a character string that identifies the type of arterial input
#' function (AIF) to be used.  Acceptable AIF models include:
#' \code{tofts.kermode}, \code{fritz.hansen}
#' @return Computed pharmacokinetic curve.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com} and Volker
#' Schmid \email{volkerschmid@@users.sourceforge.net}
#' @seealso \code{\link{dcemri.lm}}, \code{\link{dcemri.bayes}},
#' \code{\link{dcemri.spline}}
#' @references 
#' 
#' Fritz-Hansen, T., Rostrup, E., Larsson, H.B.W, Sondergaard, L.,
#' Ring, P. and Henriksen, O. (1993) Measurement of the arterial concentration
#' of Gd-DTPA using MRI: A step toward quantitative perfusion imaging,
#' \emph{Magnetic Resonance in Medicine}, \bold{36}, 225-231.
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
#' xi <- seq(5, 300, by=5)
#' img <- array(t(breast$data)[,xi], c(13,1,1,60))
#' mask <- array(TRUE, dim(img)[1:3])
#' time <- buckley$time.min[xi]
#' 
#' fit.lm <- dcemri.lm(img, time, mask, aif="fritz.hansen")
#' par.lm <- c("vp"=fit.lm$vp[3], "ktrans"=fit.lm$ktrans[3], "kep"=fit.lm$kep[3])
#' curve.lm <- kineticModel(time, par.lm)
#' plot(time, img[3,1,1,], xlab="time", ylab="contrast agent concentration")
#' lines(time, curve.lm, lwd=2, col=2)
#' 
#' fit.bayes <- dcemri.bayes(img, time, mask, aif="fritz.hansen")
#' par.bayes <- c("vp"=fit.bayes$vp[3], "ktrans"=fit.bayes$ktrans[3],
#'                "kep"=fit.bayes$kep[3])
#' curve.bayes <- kineticModel(time, par.bayes)
#' lines(time, curve.bayes, lwd=2, col=4)
#' legend("bottomright", c("Levenburg-Marquardt (extended/fritz.hansen)",
#'                         "Bayesian Estimation (extended/fritz-hansen)"),
#'        lwd=2, col=c(2,4))
#' cbind(time, img[3,,,], curve.lm, curve.bayes)[20:30,]
#' 
#' @export kineticModel
kineticModel <- function(time, par, model="extended", aif="fritz.hansen") {

  d <- dim(par["ktrans"])
  if (!is.numeric(d)) {
    d <- length(par["ktrans"])
  }
  T <- length(time)
  dd <- prod(d)
  
  if (!(is.numeric(par["kep"]))) {
    par["kep"] <- par["ktrans"] / par["ve"]
  }
  
  p <- aifParameters(aif)
  func.model <- compartmentalModel(model)
  result <- switch(model,
                   weinmann = {
                     func.model(rep(time, dd),
                                ##rep(log(par$ktrans), each=TRUE),
                                ##rep(log(par$kep), each=TRUE),
                                rep(log(par), each=TRUE),
                                p)
                     },
                   extended = {
                     func.model(rep(time, dd),
                                ##rep(log(par$vp), each=TRUE),
                                ##rep(log(par$ktrans), each=TRUE),
                                ##rep(log(par$kep), each=TRUE),
                                rep(log(par), each=TRUE),
                                p)
                   },
                   stop("Model is not currently supported."))
  result <- array(result, c(T,dd))
  result <- aperm(result, c(2:length(dim(result)),1))
  return(drop(result))
}

