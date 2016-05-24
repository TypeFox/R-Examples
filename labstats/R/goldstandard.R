#' @title Data to compare a new method against a gold standard
#'
#' @description Simulated data for one gold standard and four other
#' outcome variables.
#'
#' @details The data contains a gold standard plus four other
#' simulated variables that vary from the gold standard in different
#' ways.
#'
#' @format A data frame with 20 rows and 5 variables:
#' \describe{
#'   \item{goldstandard:}{Values for the gold standard outcome variable.}
#' 
#'   \item{perfect.cor:}{An outcome variable that is identical to the
#' gold standard.}
#' 
#'   \item{noise:}{An outcome variable created by adding noise to the
#' gold standard.}
#' 
#'   \item{location:}{An outcome variable created by adding noise plus
#' a shift in location.}
#' 
#'   \item{scale:}{An outcome variable created by adding noise plus a
#' shift in scale.}  }
#'
#' @docType data
#' @name goldstandard
NULL
