#' Rexperigen: an R interface to Experigen
#'
#' An R interface for downloading results from an Experigen server.
#' Works with the "classic" server that is currently running on
#' \code{db.experigen.org} as well, but its main advantage is that it helps
#' a lot with the new functions of the newer version of the Experigen
#' server. Most importantly, it helps with registration of experimenters,
#' registration of experiments and accessing their data.
#'
#' @section Setup functions:
#'
#' \code{\link{setExperigenServer}}, \code{\link{setExperigenCredentials}}
#'
#' @section Registration functions:
#'
#' \code{\link{registerExperiment}}, \code{\link{removeRegistration}}, \code{\link{getRegisteredExperiments}}
#'
#' @section Download functions:
#'
#' The main function is \code{\link{downloadExperiment}}. Further information
#' about the experiment can be inquired with \code{\link{getUsers}} and
#' \code{\link{getDestinations}}
#'
#' @importFrom utils URLencode read.table
#'
#' @docType package
#' @name Rexperigen
NULL

