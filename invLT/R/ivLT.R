#Christopher Barry, last updated 16/07/2015 at University of Birmingham

#Numerical inversions of the inverse Laplace Transform, with Evans & Chung's (2000) optimum contour (1st one) and the Bromwich contour

#The Evans & Chung (2000) optimum contour----

#' @title Optimum Contour
#' @description The optimum contour with polar co-ordinates (\eqn{r} as a function of \eqn{\phi}) and complex length increment with \eqn{\phi} (Evans and Chung, 2000)
#'
#' @param phi \eqn{\phi} value
#' @param m width of the contour - too small and get too close to singularities on negative \eqn{x}-axis, too large and encounter instability for large positive \eqn{x}
#' @param t standard (time) domain variable, also affects contour width
#' @details if \code{t} is set as zero, it is changed to 5 (avoids dividing by 0)
#'
#' @references Evans & Chung, 2000: Laplace transform inversions using optimal contours in the complex plane \emph{International Journal of Computer Mathematics} v73 pp531-543.
#'
#' @describeIn opC.r \eqn{r}
opC.r <- function(phi, m=1, t=5){ #optimal contour expressed in polar co-ordinates
  if(identical(t, 0)){t <- 5} #note value at t = 0 may still be wrong
  return(m*phi/(t*sin(phi)))
}

#' @describeIn opC.r \eqn{ds/d\phi}
opC.ds_dphi <- function(phi, m=1, t=5){ #derivative of contour path length with respect to phi, analytical
  #return(((m*phi/(t*sin(phi)))^2 + (1/sin(phi) - phi/(tan(phi)*sin(phi)))^2)^0.5) #s as magnitude: real
  #s as dx and dy components: complex
  if(identical(t, 0)){t <- 100} #note value at t = 0 may still be wrong; this avoids dividing by 0
  if(identical(phi, 0)){
    dx_dphi <- 0 #I have checked this by investigating small values of phi
  }else{
    dx_dphi <- (m/t)*((sin(phi) - phi*cos(phi))*cos(phi)/(sin(phi))^2 - phi)
  }
  dy_dphi <- (m/t) #independent of phi
  return(dx_dphi + dy_dphi*1i)
}

#' @title Inverse Laplace Transform
#' @description Functionals that numerically invert a Laplace Transform.
#'
#' @param L.FUN the Laplace-Transformed function
#' @param t standard (time) domain function at which to evaluate
#' @param nterms number of terms to use in the numerical inversion (odd number safest for \code{iv.opC}, even for \code{iv.opChalf})
#' @param m see \code{\link{opC.r}} documentation
#' @param fail.val value to return in event of failure to converge
#' @param gamma the Bromwich contour is a straight line and intersects the real axis at \eqn{\gamma}
#'
#' @details Optimum contour based on:
#'
#' Evans & Chung, 2000: Laplace transform inversions using optimal contours in the complex plane \emph{International Journal of Computer Mathematics} v73 pp531-543.
#'
#' @examples
#' tvals <- seq(-pi/2, pi/2, length.out = 7)
#' sinvals <- vapply(tvals, iv.opC, complex(1), L.FUN = L.sin)
#' plot(tvals, Re(sinvals), type = "l")
#'
#' @describeIn iv.opC inversion using the full optimum contour
iv.opC <- function(L.FUN, t, nterms = 31L, m=1, fail.val = NA){
  FUN <- match.fun(L.FUN) #ensures that FUN is read as a function
  if(identical(t, 0)){t <- 10^-50} #accurate for very low values of t, but not exactly 0
  for(attempt in 1:10){
    dphi <- 2*pi/nterms
    phi <- -pi + (1:nterms - 1/2)*dphi
    z <- opC.r(phi, m, t)*cos(phi) + opC.r(phi, m, t)*sin(phi)*1i
    L <- vapply(z, L.FUN, 1i)
    if(any(!is.finite(L))){
      if(identical(attempt, 10L)){cat("Laplace Transform inversion failed after 10 attempts.\n"); return(fail.val)}
      m <- m*2
      nterms <- round(nterms*1.378, 0) #random irregular number to avoid sampling same point again
      next
    }
    em <- exp(z*t)
    ds_dphi <- vapply(phi, opC.ds_dphi, 1i, m, t)
    break
  }
  return(sum(em*L*ds_dphi)*dphi/(2i*pi))
}

#' @describeIn iv.opC for functions which are symmetric about the real axis, it is sufficient to use half the optimum contour and half the number of subdivisions (\code{nterms})
iv.opChalf <- function(L.FUN, t, nterms = 16L, m=1, fail.val = NA){
  FUN <- match.fun(L.FUN) #ensures that FUN is read as a function
  if(identical(t, 0)){t <- 10^-50} #accurate for very low values of t, but not exactly 0
  for(attempt in 1:10){
    dphi <- pi/nterms
    phi <- (1:nterms - 1/2)*dphi
    z <- opC.r(phi, m, t)*cos(phi) + opC.r(phi, m, t)*sin(phi)*1i
    L <- vapply(z, L.FUN, 1i)
    if(any(!is.finite(L))){
      if(identical(attempt, 4L)){cat("Laplace Transform inversion failed after", attempt, "attempts, with t =", signif(t, 4), ". ", fail.val, "returned.\n"); return(fail.val)} # have to stop somewhere - return whatever value was specified in case of failure: NA is used for safety
      m <- m*2
      nterms <- round(nterms*1.378, 0) #random irregular number to avoid sampling same point again
      next
    }
    em <- exp(z*t)
    ds_dphi <- vapply(phi, opC.ds_dphi, 1i, m, t)
    break
  }
  return(sum(Im(em*L*ds_dphi))*dphi/pi)
}

#the Bromwich contour----

#' @title Bromwich Contour
#' @description The Bromwich contour with polar co-ordinates (\eqn{r} as a function of \eqn{\phi})
#'
#' @param phi \eqn{\phi} value
#' @inheritParams iv.opC
#'
#' @describeIn BrC.r \eqn{r}
BrC.r <- function(phi, gamma = 1){return(gamma/cos(phi))} #Bromwich contour

#' @describeIn BrC.r \eqn{ds/d\phi}
BrC.ds_dphi <- function(phi, gamma = 1){
  return(1i*(BrC.r(phi, gamma)^2 + (gamma*tan(phi)/cos(phi))^2)^0.5) #purely imaginary
}

#' @describeIn iv.opC inversion using the Bromwich contour (the definition, but very unstable for numerical evaluation - not recommended)
iv.BrC <- function(L.FUN, t, nterms = 1000L, gamma = 1){
  FUN <- match.fun(L.FUN) #ensures that FUN is read as a function
  if(t == 0){t <- 10^-50}
  tot <- 0 #initialise
  dphi <- pi/nterms
  for(n in 1:nterms){
    phi <- -pi/2 + (n - 1/2)*dphi
    ds <- dphi*BrC.ds_dphi(phi, gamma)
    x <- x.rphi(BrC.r(phi, gamma), phi)
    y <- y.rphi(BrC.r(phi, gamma), phi)
    tot <- tot + exp((x + y*1i)*t)*FUN(x + y*1i)*ds
  }
  return(tot/(2*pi*1i))
}

#Tools----

#' @title Plot Laplace Transform inversion
#' @description Plots the results of a Laplace Transform inversion at multiple time values.
#'
#' @inheritParams iv.opC
#' @param METHOD inversion algorithm to use (iv.opC, iv.opChalf or iv.BrC)
#' @param tPts time points at which to plot
#' @param ... graphical parameters for \code{\link[graphics]{plot}}
#'
#' @details This function is useful for investigating the performance of a Laplace Transform inversion over a range of time values.  Use for example with the LT functions provided in with this package (invLT).
#'
#' @examples
#' ivLT.plot(L.tsq, iv.opC, nterms = 31L)
#' ivLT.plot(L.tsq, iv.opC, nterms = 1000L)
#' ivLT.plot(L.tsq, iv.opChalf, nterms = 16L)
#' ivLT.plot(L.tsq, iv.opChalf, nterms = 1000L)
#' ivLT.plot(L.tsq, iv.BrC, nterms = 31L)
#' ivLT.plot(L.tsq, iv.BrC, nterms = 1000L)
ivLT.plot <- function(L.FUN, METHOD = iv.opC, tPts = seq(-2,5,.1), nterms = 100, ...){
  FUN <- match.fun(L.FUN)
  METHOD <- match.fun(METHOD)
  iv.LT <- rep(0, length(tPts))
  for(tn in 1:length(tPts)){iv.LT[tn] <- METHOD(FUN, tPts[tn], nterms = nterms)}
  plot(tPts, Re(iv.LT), type = "l", xlab = "t", ylab = "inv. LT", ...)
}

#' @title Laplace Transforms
#' @description Laplace Transforms of common functions.  Useful for testing out LT inversion functions and whether sufficient precision is being used.
#'
#' @param p Laplace domain variable (commonly called \eqn{s} elsewhere)
#'
#' @describeIn L.t LT of \eqn{t}
L.t <- function(p){1/p^2}# LT of t^1

#' @describeIn L.t LT of \eqn{t^2}
L.tsq <- function(p){2/p^3}# LT of t^2

#' @describeIn L.t LT of \eqn{e^(-t)}
L.exp <- function(p){1/(p+1)}# LT of exp(-t)

#' @describeIn L.t LT of cos\eqn{(t)}
L.cos <- function(p){p/(p^2 + 1)}# LT of cos(t)\n

#' @describeIn L.t LT of sin\eqn{(t)}
L.sin <- function(p){1/(p^2 + 1)}# LT of sin(t)\n

#' @describeIn L.t LT of Heaviside unit function stepping at 1: (if \eqn{p < 1} 0 else 1)
L.H <- function(p){exp(-p)/p}# LT of H(t - 1), the Heaviside unit function stepping at t = 1\n

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/Scripts/R",
    devtools.install.args = "",
    devtools.name = "invLT",
    devtools.desc.author = '"Christopher Barry <cjb309@bham.ac.uk> [Barry, 2015]"',
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])

  invisible()
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage("Numerical Laplace Transform inversion functions successfully sourced.\n",
      "Optimum contour integration ref. Evans & Chung, 2000: Laplace transform inversions using optimal contours in the complex plane; International Journal of Computer Mathematics v73 pp531-543.:\n",
      "iv.opC(L.FUN, t, nterms = 30L, m=1, fail.val = NA){\n",
      "iv.opC2(L.FUN, t, nterms = 30L, m=1, fail.val = NA){\n",
      "L.FUN: Laplace Transformed function to be inverted\n",
      "t: untransformed domain co-ordinate value (typically time) at which to evaluate\n",
      "nterms: number of terms with which to evaluate the inversion\n",
      "m: \"width\" of the contour of integration - smaller values avoid the large oscillations at the right hand side of the Argand diagram, but are more likely to interact with singlurities on the real axis\n",
      "for further details on m, see E&C 2000, or simply experiment\n",
      "fail.val: after 10 attempts, increasing m and nterms, if the inversion still produces non-finite results, the function gives up and returns this value\n\n",
      "iv.opC2 only uses half of the contour and therefore nterms can be half the size for the same accuracy: valid if Re(L.FUN) is symmetric about the real axis\n\n",
      "Bromwich contour integration:\n",
      "iv.BrC(L.FUN, t, nterms=1000, gamma=1)\n",
      "gamma: distance of contour to right of imaginary axis\n\n",
      "The Bromwich Contour is very unstable; it is included mainly for comparison.  The definition of the inverse Laplace Transform is usually expressed with the Bromwich Contour, even though its numerical implementation is impractical.\n\n",
      "Tools:\n",
      "ivLT.plot(L.FUN, METHOD = iv.opC, tPts = seq(-2,5,.1), nterms = 100, ...)\n",
      "METHOD: iv.opC, iv.opC2 or iv.BrC\n",
      "tPts: points at which to plot time\n",
      "...: parameters to pass to plot\n",
      "some transformed functions for which there are analytical inverses:\n",
      "L.t <- function(p){1/p^2}: LT of t^1\n",
      "L.tsq <- function(p){2/p^3}: LT of t^2\n",
      "L.exp <- function(p){1/(p+1)}: LT of exp(-t)\n",
      "L.cos <- function(p){p/(p^2 + 1)}: LT of cos(t)\n",
      "L.sin <- function(p){1/(p^2 + 1)}: LT of sin(t)\n",
      "L.H <- function(p){exp(-p)/p}: LT of H(t - 1), the Heaviside unit function stepping at t = 1\n")
}
