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

#' @title Calculate the Asymptotic Covariance Matrix
#' @description Places the Asymptotic Covariance Matrix in print form.
#' @param signal.modwt A \code{modwt} object that contains the modwt decomposition.
#' @param signal.wvar  A \code{wvar} object that contains the wavelet variance.
#' @param compute.v    A \code{string} that indicates the type of covariance matrix to compute. Supports: "diag"
#' @return A \code{list} with the structure:
#' \describe{
#'   \item{"V"}{Covariance Matrix}
#'   \item{"V.r"}{Covariance Matrix}
#'   \item{"nlevels"}{Level of decomposition J}
#'   \item{"compute.v"}{Type of Covariance Matrix}
#'   \item{"robust"}{Robust active}
#'   \item{"eff"}{Efficiency level for Robust}
#'   \item{"scales"}{Tau scales (2^(1:J))}
#'   \item{"wv.empir"}{Empirical Wavelet Variance}
#' }
#' @author JJB
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' decomp = modwt(x)
#' wv = wvar(x)
#' out = wvcov(decomp, wv, compute.v="diag")
#' 
#' # Robust
#' decomp = modwt(x)
#' wv = wvar(x, robust = TRUE)
#' out = wvcov(decomp, wv, compute.v="diag")
wvcov = function(signal.modwt, signal.wvar, compute.v="diag"){
  
  if(!is(signal.modwt,"modwt")){
    stop("Need to supply a gmwm_modwt object as the first parameter.")
  }
  
  if(!is.wvar(signal.wvar)){
    stop("Need to supply a wvar object as the second parameter.")
  }
  
  out = .Call('gmwm_compute_cov_cpp', PACKAGE = 'gmwm', signal.modwt, attr(signal.modwt,"J"), compute.v, signal.wvar$robust, signal.wvar$eff)
  out = structure(list(V=out[[1]],
                       V.robust=out[[2]], 
                       nlevels=attr(signal.modwt,"J"), 
                       compute.v = compute.v, 
                       robust = signal.wvar$robust, 
                       eff = signal.wvar$eff, 
                       scales = signal.wvar$scales,
                       wv.empir = signal.wvar$variance,
                       ci_low = signal.wvar$ci_low,
                       ci_high = signal.wvar$ci_high), class = "wvcov")
  invisible(out)
}



#' @title Print Asymptotic Covariance Matrix
#' @description Places the Asymptotic Covariance Matrix in print form.
#' @method print wvcov
#' @export
#' @keywords internal
#' @param x    A \code{wvcov} object
#' @param ...  Further arguments passed to or from other methods
#' @return Prints the modwt matrix decomposition
#' @author JJB
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' decomp = modwt(x)
#' wv = wvar(x)
#' print(wvcov(decomp,wv,compute.v="diag"))
print.wvcov = function(x, ...){
  print(x$V)
}

#' @title Summary Wavelet Covariance Matrix
#' @description Prints the Wavelet Covariance Matrix
#' @method summary wvcov
#' @export
#' @keywords internal
#' @param object A \code{wvcov} object
#' @param ...    Additional arguments affecting the summary produced.
#' @return Prints the modwt matrix decomposition
#' @author JJB
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' decomp = modwt(x)
#' wv = wvar(x)
#' summary(wvcov(decomp,wv,compute.v="diag"))
summary.wvcov=function(object, ...){
  
  name = if(object$robust){
    "robust" 
  }else{
    "classical"
  }
  cat("The asymptotic ", object$compute.v, " using the ",name, " method.\n",sep="")
  if(object$robust){
    cat("The robust form was created using efficiency=",object$eff,"\n",sep="")
  }  
  
  print(object)
}