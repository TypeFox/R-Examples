#' Interface between R and vcglib libraries for mesh operations
#' 
#' Provides meshing functionality from vcglib (meshlab) for R. E.g. mesh
#' smoothing, mesh decimation, closest point search.
#' 
#' \tabular{ll}{
#' Package: \tab Rvcg\cr
#' Type: \tab Package\cr
#' Version: \tab 0.13.1.2\cr
#' Date: \tab 2016-01-21\cr
#' License: \tab GPL\cr
#' LazyLoad: \tab yes\cr }
#' 
#' @name Rvcg-package
#' @aliases Rvcg-package Rvcg
#' @docType package
#' @author Stefan Schlager
#' 
#' Maintainer: Stefan Schlager <zarquon42@@gmail.com>
#' @references To be announced
#' @keywords package
#' @import grDevices stats utils
#' @importFrom Rcpp evalCpp 
#' @useDynLib Rvcg
NULL

#' Example mesh and landmarks
#'
#' A triangular mesh representing a human face - called by data(humface)
#' 
#' @name humface
#' @aliases humface humface.lm
#' @docType data
#' @format \code{humface}: triangular mesh representing a human face.
#'
#' \code{humface.lm}: landmarks on mesh 'humface'- called by data(humface)
#'
#' @keywords datasets

NULL


#' dummyhead - dummy head and landmarks
#'
#' A triangular mesh representing a dummyhead - called by data(dummyhead)
#' 
#' @name dummyhead
#' @aliases dummyhead.mesh dummyhead.lm
#' @docType data
#' @format \code{dummyhead.mesh}: triangular mesh representing a dummyhead.
#'
#' \code{dummyhead.lm}: landmarks on mesh 'dummyhead'
#' @keywords datasets
#' 
NULL
