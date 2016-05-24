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

#' @title Calculate the Allan Variance
#' @description Computes the Allan Variance
#' @param x     A \code{vec} containing the time series under observation.
#' @param type  A \code{string} containing either \code{"mo"} for Maximal Overlap or \code{"to"} for Tau Overlap
#' @return Allan variance fixed
#' @author JJB
#' @return av A \code{list} that contains:
#' \itemize{
#'  \item{"clusters"}{The size of the cluster}
#'  \item{"allan"}{The Allan variance}
#'  \item{"errors"}{The error associated with the variance estimation.}
#' }
#' @details 
#' The decomposition and the amount of time it takes to perform it depends on whether you are using
#' the Tau Overlap or the Maximal Overlap.
#' 
#' Maximal Overlap Allan Variance
#' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
#' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
#' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n|n< floor(log2(N))}}
#' Then, \eqn{M = N - 2n} samples exist. 
#' The Maximal-overlap estimator is given by:
#' \eqn{\frac{1}{{2\left( {N - 2k + 1} \right)}}\sum\limits_{t = 2k}^N {{{\left[ {{{\bar Y}_t}\left( k \right) - {{\bar Y}_{t - k}}\left( k \right)} \right]}^2}} }
#' 
#' where \eqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }.
#' 
#' Tau-Overlap Allan Variance
#' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
#' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
#' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n|n< floor(log2(N))}}
#' Then, a sampling of \eqn{m = \left\lfloor {\frac{{N - 1}}{n}} \right\rfloor  - 1} samples exist. 
#' The tau-overlap estimator is given by:
#' 
#' where \eqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }.
#' 
#' @author JJB
#' @references Long-Memory Processes, the Allan Variance and Wavelets, D. B. Percival and P. Guttorp
#' @examples
#' set.seed(999)
#' # Simulate white noise (P 1) with sigma^2 = 4
#' N = 100000
#' white.noise = rnorm(N, 0, 2)
#' #plot(white.noise,ylab="Simulated white noise process",xlab="Time",type="o")
#' #Simulate random walk (P 4)
#' random.walk = cumsum(0.1*rnorm(N, 0, 2))
#' combined.ts = white.noise+random.walk
#' av_mat = avar_mo_cpp(combined.ts)
avar = function(x, type = "mo") {
  x = as.vector(x)
  
  if(type == "mo"){
    av = .Call('gmwm_avar_mo_cpp', PACKAGE = 'gmwm', x)
  }else{
    av = .Call('gmwm_avar_to_cpp', PACKAGE = 'gmwm', x)
  }
  
  av = list(clusters = av[,1], allan=av[,2], errors=av[,3])
  av$adev = sqrt(av$allan)
  av$lci = av$adev - av$errors*av$adev
  av$uci = av$adev + av$errors*av$adev
  class(av) = "avar"
  av
}

#' @title Prints Allan Variance
#' @description Displays the allan variance information
#' @method print avar
#' @export
#' @param x   A \code{avar} object.
#' @param ... Arguments to be passed to methods
#' @author JJB
#' @return console output
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = avar(x)
#' print( out )
print.avar = function(x, ...) {
  cat("\n Clusters: \n")
  print(x$clusters, digits=5)
  cat("\n Allan Variances: \n")
  print(x$allan, digits=5)
  cat("\n Errors: \n")
  print(x$errors, digits=5)
}

#' @title Plot Allan Variance
#' @description Displays a plot containing the allan variance
#' @method plot avar
#' @export
#' @param x   A \code{avar} object.
#' @param ... Arguments to be passed to methods
#' @author JJB
#' @return ggplot2 graph
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = avar(x)
#' plot( out )
plot.avar = function(x, ...){
  plot(x$clusters, x$adev,log="xy",
       xlab=expression("Scale " ~ tau),
       ylab=expression("Allan Deviation " ~ phi[tau]),
       main=expression(log(tau) ~ " vs. " ~ log(phi[tau]))
  )
  lines(x$clusters, x$adev, type="l")
  lines(x$clusters, x$lci, type="l", col="grey", lty=2)
  lines(x$clusters, x$uci, type="l", col="grey", lty=2)
}


#' @title Summary Allan Variance
#' @description Displays the summary table of allan variance
#' @method summary avar
#' @export
#' @param object A \code{avar} object.
#' @param ...    Additional arguments affecting the summary produced.
#' @author JJB
#' @return Summary table
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = avar(x)
#' summary( out )
summary.avar = function(object, ...) {
  out_matrix = matrix(0, nrow = length(object$clusters), ncol = 6)
  colnames(out_matrix) = c("Time", "AVAR", "ADEV", "Lower CI", "Upper CI", "Error")
  out_matrix[,"Time"] = object$clusters
  out_matrix[,"AVAR"] = object$allan
  out_matrix[,"ADEV"] = object$adev
  out_matrix[,"Lower CI"] = object$lci
  out_matrix[,"Upper CI"] = object$uci
  out_matrix[,"Error"] = object$errors
  
  class(out_matrix) = "summary.avar"
  out_matrix
}