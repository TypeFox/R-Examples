#' Smooth Data in R
#' 
#' @description \code{smoother} Package for the Smoothing of Numerical Data
#' 
#' @details \code{smoother} is presently limited to a port of the Matlab 'Gaussian Window' Function, 
#' as well as a limited number of moving averages (\code{sma, ema, dema} and \code{'wma'}). Code for the gaussian window 
#' function has been written locally within this package, however, the moving averages are called from the \link{TTR} package 
#' (\url{http://cran.r-project.org/web/packages/TTR/index.html}) and are included as a matter of convenience.
#' 
#' For further information (and examples) with regards to utilizing the primary helper function, 
#' please refer to the \link{smth} function help file
#' 
#' @references The Gaussian Smoothing component of the \code{smoother} package has been loosley adapted from the 
#' following works: \url{http://goo.gl/NK79bJ}.
#' 
#' @name smoother
#' @rdname smoother
#' @aliases smoother
NULL
