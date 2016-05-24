# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.

#' @title Brickwall functionality for MO/DWT
#' @description 
#' Removes boundary coefficients
#' @param signal.decomp  A \code{modwt} or \code{dwt} object that has not yet been brick walled
#' @return A \code{modwt} or \code{dwt} object that has been brick walled
#' @author JJB
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' o = modwt(x, bw = FALSE)
#' brickwall(o)
#' 
#' x = rnorm(2^8)
#' j = dwt(x, bw = FALSE)
#' brickwall(j)
brickwall = function(signal.decomp){
  if(!(is(signal.decomp,"modwt") || is(signal.decomp,"dwt"))){
    stop("`signal.decomp` must be from the either the `modwt()` or `dwt()` function.")
  }
  if(attr(signal.decomp,"brick.wall")){
    stop("The decomposition has already been decomposed.")
  }
  
  obj = .Call('gmwm_brick_wall', PACKAGE = 'gmwm', signal.decomp, select_filter("haar"), class(signal.decomp))
  
  mostattributes(obj) = attributes(signal.decomp)
  attr(signal.decomp,"brick.wall") = T
  
  obj
}



#' @title Maximum Overlap Discrete Wavelet Transform
#' @description 
#' Calculation of the coefficients for the discrete wavelet transformation
#' @param x        A \code{vector} with dimensions N x 1. 
#' @param nlevels  A \code{integer} indicating the \eqn{J} levels of decomposition.
#' @param filter   A \code{string} indicating the filter name
#' @param boundary A \code{string} indicating whether the filter is: \code{"periodic"} or \code{"reflection"}.
#' @param bw       A \code{boolean} indicating whether to remove (TRUE) or keep (FALSE) boundary wavelet coefficients
#' @return y       A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
#' @details
#' Performs a level \eqn{J} decomposition of the time series using the pyramid algorithm.
#' The default \eqn{J} is determined by \eqn{floor\left(log_2 \left(length\left(x\right)\right)\right)}{floor(log2(length(x)))}
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' modwt(x)
modwt = function(x, nlevels = floor(log2(length(x))), filter = "haar", boundary="periodic", bw = TRUE) {
  out = .Call('gmwm_modwt_cpp', PACKAGE = 'gmwm', x, filter_name = filter, nlevels, boundary = boundary, brickwall = bw)
  names(out) = paste0("S",1:nlevels)
  mostattributes(out) = list(J=nlevels, filter = filter, boundary = boundary, brick.wall = bw, class="modwt")
  out
}

#' @title Print Maximum Overlap Discrete Wavelet Transform
#' @description
#' Prints the results of the modwt list
#' @method print modwt
#' @export
#' @param x A \code{modwt} object
#' @param ... further arguments passed to or from other methods.
#' @return Prints the modwt decomposition
#' @author JJB
#' @keywords internal
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' print(modwt(x))
print.modwt=function(x, ...){
  x
}

#' @title Summary Maximum Overlap Discrete Wavelet Transform
#' @description Unlists MODWT object and places it in matrix form
#' @method summary modwt
#' @export
#' @keywords internal
#' @param object A \code{modwt} object
#' @param ... additional arguments affecting the summary produced.
#' @return Prints the modwt matrix decomposition
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' summary(modwt(x))
summary.modwt=function(object, ...){
  cat("Results of the MODWT containing ",attr(object,"J")," scales\n")
  cat("These values are", if(!attr(object,"brick.wall")){" >NOT<"}," brick walled\n")
  print.modwt(object)
}