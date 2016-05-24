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

#' @title Discrete Wavelet Transform
#' @description 
#' Calculation of the coefficients for the discrete wavelet transformation
#' @inheritParams modwt
#' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
#' @details
#' Performs a level \eqn{J} decomposition of the time series using the pyramid algorithm.
#' The default \eqn{J} is determined by \eqn{floor\left(log_2 \left(length\left(x\right)\right)\right)}{floor(log2(length(x)))}
#' @author JJB
#' @examples
#' set.seed(999)
#' x = rnorm(2^8)
#' dwt(x)
dwt = function(x, nlevels = floor(log2(length(x))), filter = "haar", boundary="periodic", bw = TRUE) {
  
  if(is.vector(x) && length(x) %% nlevels != 0){
    warning("The data has been truncated so that it is divisible by `nlevels` (e.g. 2^*)")
    x = x[1:2^nlevels]
  }else if(is.matrix(x) || is.data.frame(x)){
    if(ncol(x) != 1){
      stop("Only one column is allowed to be decomposed at a time.")
    }
    
    if(nrow(x) %% nlevels !=0){
      warning("The data has been truncated so that it is divisible by `nlevels` (e.g. 2^*)")
      idx = 1:2^nlevels
      x[idx,1] = x[idx,1]
    }
  }
  out = .Call('gmwm_dwt_cpp', PACKAGE = 'gmwm', x, filter_name = filter, nlevels, boundary = boundary, brickwall = bw)
  names(out) = paste0("S",1:nlevels)
  mostattributes(out) = list(J=nlevels, filter = filter, boundary = boundary, brick.wall = bw, class="dwt")
  out
}

#' @title Print Discrete Wavelet Transform
#' @description
#' Prints the results of the modwt list
#' @method print dwt
#' @export
#' @param x A \code{dwt} object
#' @param ... further arguments passed to or from other methods.
#' @return Prints the dwt decomposition
#' @author JJB
#' @keywords internal
#' @examples
#' set.seed(999)
#' x = rnorm(2^8)
#' print(dwt(x))
print.dwt=function(x, ...){
  x
}

#' @title Summary Discrete Wavelet Transform
#' @description Unlists DWT object and places it in matrix form
#' @method summary dwt
#' @export
#' @keywords internal
#' @param object A \code{dwt} object
#' @param ... additional arguments affecting the summary produced.
#' @return Prints the dwt matrix decomposition
#' @author JJB
#' @examples
#' set.seed(999)
#' x = rnorm(2^8)
#' summary(dwt(x))
summary.dwt=function(object, ...){
  cat("Results of the DWT containing ",attr(object,"J")," scales\n")
  cat("These values are", if(!attr(object,"brick.wall")){" >NOT<"}," brick walled\n")
  print.dwt(object)
}