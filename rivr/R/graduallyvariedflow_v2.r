#' Steady and Unsteady Open-Channel Flow Computation
#'
#' This package is designed as an educational tool for students and 
#' instructors of undergraduate courses in open channel hydraulics. 
#' Functions are provided for computing normal and critical depths, 
#' steady (e.g. backwater curves) and unsteady (flood wave routing) 
#' flow computations for prismatic trapezoidal channels. See the vignettes
#' to get started.
#' @name rivr-package
#' @aliases rivr
#' @docType package
#' @useDynLib rivr
#' @importFrom Rcpp evalCpp
NULL

get_profile = function(So, n, Q, g, Cm, B, SS, y0){
  yc = critical_depth(Q, y0, g, B, SS)
  yn = normal_depth(So, n, Q, y0, Cm, B, SS)
  if(So < 0) # adverse slope
    if (y0 > yc)
      return("A2")
    else
      return("A3")
  else if (So == 0) # horizontal slope
    if(y0 > yc)
      return("H2")
    else
      return("H3")
  else if (yn > yc) # Mild slope
    if(y0 > yn)
      return("M1")
    else if (y0 > yc)
      return("M2")
    else
      return("M3")
  else if (yc > yn) # steep slope
    if (y0 > yc)
      return("S1")
    else if (yn > y0)
      return("S3")
    else
      return("S2")    
  else # critical profile
    if (y0 > yc)
      return("C1")
    else
      return("C3")    
}

check_profile = function(p){
  if(any(p == c("S1", "M3", "H3", "A3"))){
    stop(p, " profile cannot be computed (rapidly-varied flow)")  
  } else if(any(p == c("A2", "H2", "M1", "M2", "C1"))){
    message(p, " profile specified. Computing upstream profile")
    function(x) -abs(x)
  } else {
    message(p, " profile specified. Computing downstream profile")
    function(x) abs(x)
  } 
}

#' @title Gradually-varied flow profiles
#' @description Compute the gradually-varied flow profile of a prismatic channel.
#' @param So Channel slope [\eqn{L L^{-1}}].
#' @param n Manning's roughness coefficient.
#' @param Q Flow rate [\eqn{L^3 T^{-1}}].
#' @param y0 The water depth at the control section [\eqn{L}].
#' @param Cm Unit conversion coefficient for Manning's equation. For SI units, Cm = 1.
#' @param g Gravitational acceleration [\eqn{L T^{-2}}].
#' @param B Channel bottom width [\eqn{L}].
#' @param SS Channel sideslope [\eqn{L L^{-1}}].
#' @param z0 Elevation reference datum at control section [\eqn{L}]. Default is 0.
#' @param x0 Distance reference at control section [\eqn{L}]. Default is 0.
#' @param stepdist The spatial interval used in the Standard step method [\eqn{L}].
#' @param totaldist The total distance upstream (or downstream) to compute the profile [\eqn{L}].
#' @return data.frame with columns:
#'   \item{x}{Along-channel distance.}
#'   \item{z}{Elevation.}
#'   \item{y}{Flow depth.}
#'   \item{v}{Flow velocity.}
#'   \item{A}{Flow area.}
#'   \item{Sf}{Friction slope.}
#'   \item{E}{Total energy.}
#'   \item{Fr}{Froude Number.}
#' @details Computes the longitudinal water surface profile of a prismatic 
#'   channel using the standard step method by solving the non-linear ODE 
#'   \deqn{\frac{dy}{dx} = \frac{S_0 - S_f}{1 - Fr^2}} The standard-step 
#'   method operates by stepping along the channel by a constant distance 
#'   interval, starting from a cross-section where the flow depth is known 
#'   (the control section). The flow depth is computed at the adjacent 
#'   cross-section (target section). The computed value at the target is then 
#'   used as the basis for computing flow depth at the next cross-section, i.e. 
#'   the previous target section becomes the new control section for each step. 
#'   A Newton-Raphson scheme is used each step to compute the flow depth and 
#'   friction slope. Technically, the average friction slope of the control and
#'   target section is used to compute the flow depth at the target section.
#' @examples
#' # example M1 profile
#' compute_profile(0.001, 0.045, 250, 2.7, 1.486, 32.2, 100, 0, stepdist = 10, totaldist = 3000)
#' # example M2 profile
#' compute_profile(0.001, 0.045, 250, 0.64, 1.486, 32.2, 100, 0, stepdist = 10, totaldist = 3000)
#' # example S2 profile
#' compute_profile(0.005, 0.01, 250, 2.65, 1.486, 32.2, 10, 0, stepdist = 10, totaldist = 2000)
#' # example S3 profile
#' compute_profile(0.005, 0.01, 250, 0.5, 1.486, 32.2, 10, 0, stepdist = 10, totaldist = 2000)
#' @export
compute_profile = function(So, n, Q, y0, Cm, g, B, SS, z0=0, x0=0, stepdist, totaldist){
  # determine profile type
  proftype = get_profile(So, n, Q, g, Cm, B, SS, y0)
  stepsize = check_profile(proftype)(stepdist)
  res = as.data.frame(loop_step(So, n, Q, Cm, g, y0, B, SS, z0, x0, stepsize, totaldist))
  names(res) = c("x", "z", "y", "v", "A", "Sf", "E", "Fr")  
  retobj = res
  class(retobj) = c("rivr", "data.frame")  
  attr(retobj, "call") = match.call()
  attr(retobj, "simtype") = "gvf"
  attr(retobj, "proftype") = proftype
  attr(retobj, "modspec") = list(delta.x = stepdist, channel.length = totaldist,
    normal.depth = normal_depth(So, n, Q, y0, Cm, B, SS),
    critical.depth = critical_depth(Q, y0, g, B, SS))
  attr(retobj, "channel.geometry") = list(So = So, n = n, B = B, SS = SS)
  return(retobj)
}
