#Christopher Barry, started on 05/02/2015 at University of Birmingham
#last updated 08/04/2015

#co-ordinate conversion from cartesian to polar, and back

#' Conversion from cartesian to polar co-ordinates
#' @title Cartesian to Polar
#' @describeIn r.xy Returns polar co-ordinate r from cartesian co-ordinates x and y.
#'
#' @param x x co-ordinate
#' @param y y co-ordinate
#' @return r or phi respectively from x and y
r.xy <- function(x, y){
  return((x^2 + y^2)^0.5)
}

#' @describeIn r.xy Returns polar co-ordinate phi (anti-clockwise rotation from positive x-axis) from cartesian co-ordinates x and y.
phi.xy <- function(x, y){
  return(atan2(y, x)) #note use of atan2 for the sake of points left of y axis
}

#' Conversion from polar to cartesian co-ordinates
#' @title Polar to Cartesian
#' @describeIn x.rphi Returns cartesian co-ordinate x from polar co-ordinates r and phi.
#'
#' @param r distance from origin
#' @param phi anti-clockwise rotation from positive x-axis
#' @return x or t respectively from r and phi
x.rphi <- function(r, phi){
  return(r*cos(phi))
}

#' @describeIn x.rphi Returns cartesian co-ordinate x from polar co-ordinates r and phi.
y.rphi <- function(r, phi){
  return(r*sin(phi))
}
