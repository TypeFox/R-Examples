#' @import methods
NULL

#' Class "jmcmMod" of Fitted Joint Mean-Covariance Models.
#'
#' @slot call the matched call
#' @slot opt the optimization result returned by optimizeJmcm
#' @slot args arguments m, Y, X, Z, W, time
#' @slot triple an integer vector of length three containing the degrees of the
#' three polynomial functions for the mean structure, the log innovation
#' -variances and the autoregressive or moving average coefficients when
#' 'mcd' or 'acd' is specified for cov.method. It refers to the mean structure,
#' variances and angles when 'hpc' is specified for cov.method.
#' @slot devcomp the deviance components list
#'
#' @exportClass jmcmMod
setClass("jmcmMod",
  representation(
    call = "call",
    opt = "list",
    args = "list",
    triple = "numeric",
    devcomp = "list"
    ))

#' Class MCD
#'
#' Class \code{MCD} defines a joint mean covariance model based on Modified
#' Cholesky decomposition (MCD).
#'
#' @name MCD-class
#'
#' @exportClass MCD
setClass("MCD", representation(pointer="externalptr"))

MCD_method <- function(name) {
  paste("MCD", name, sep = "__")
}

#' Extract parts of MCD.
#'
#' @param x an object of MCD class
#' @param name member function of MCD class
setMethod("$", "MCD", function(x, name) {
  mcdfun <- MCD_method(name)
  function(...) .Call(mcdfun, x@pointer, ...)
})

setMethod("initialize", "MCD", function(.Object, ...) {
  .Object@pointer <- .Call(MCD_method("new"), ...)
  .Object
})

#' Class ACD
#'
#' Class \code{ACD} defines a joint mean covariance model based on Alternative
#' Cholesky decomposition (ACD).
#'
#' @name ACD-class
#'
#' @exportClass ACD
setClass("ACD", representation(pointer="externalptr"))

ACD_method <- function(name) {
  paste("ACD", name, sep = "__")
}

#' Extract parts of ACD.
#'
#' @param x an object of ACD class
#' @param name member function of ACD class
setMethod("$", "ACD", function(x, name) {
  acdfun <- ACD_method(name)
  function(...) .Call(acdfun, x@pointer, ...)
})

setMethod("initialize", "ACD", function(.Object, ...) {
  .Object@pointer <- .Call(ACD_method("new"), ...)
  .Object
})

#' Class HPC
#'
#' Class \code{HPC} defines a joint mean covariance model based on Hyperpherical
#'  parameterization of Cholesky factor (HPC).
#'
#' @name HPC-class
#'
#' @exportClass HPC
setClass("HPC", representation(pointer="externalptr"))

HPC_method <- function(name) {
  paste("HPC", name, sep = "__")
}

#' Extract parts of HPC.
#'
#' @param x an object of HPC class
#' @param name member function of HPC class
setMethod("$", "HPC", function(x, name) {
  hpcfun <- HPC_method(name)
  function(...) .Call(hpcfun, x@pointer, ...)
})

setMethod("initialize", "HPC", function(.Object, ...) {
  .Object@pointer <- .Call(HPC_method("new"), ...)
  .Object
})
