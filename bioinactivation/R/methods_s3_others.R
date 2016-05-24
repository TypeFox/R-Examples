
## "is." methods

#'
#' Test of SimulInactivation object
#'
#' Tests if an object is of class \code{SimulInactivation}.
#'
#' @param x object to be checked.
#'
#' @return A logic specifying whether \code{x} is of class
#'         \code{SimulInactivation}
#'
#' @export
#'
is.SimulInactivation <- function(x) inherits(x, "SimulInactivation")

#'
#' Test of IsoFitInactivation object
#'
#' Tests if an object is of class \code{IsoFitInactivation}.
#'
#' @param x object to be checked.
#'
#' @return A logic specifying whether \code{x} is of class
#'         \code{IsoFitInactivation}
#'
#' @export
#'
is.IsoFitInactivation <- function(x) inherits(x, "IsoFitInactivation")

#'
#' Test of FitInactivation object
#'
#' Tests if an object is of class \code{FitInactivation}.
#'
#' @param x object to be checked.
#'
#' @return A logic specifying whether \code{x} is of class
#'         \code{FitInactivation}
#'
#' @export
#'
is.FitInactivation <- function(x) inherits(x, "FitInactivation")

#'
#' Test of FitInactivationMCMC object
#'
#' Tests if an object is of class \code{FitInactivationMCMC}.
#'
#' @param x object to be checked.
#'
#' @return A logic specifying whether \code{x} is of class
#'         \code{FitInactivationMCMC}
#'
#' @export
#'
is.FitInactivationMCMC <- function(x) inherits(x, "FitInactivationMCMC")

#'
#' Test of PredInactivationMCMC object
#'
#' Tests if an object is of class \code{PredInactivationMCMC}.
#'
#' @param x object to be checked.
#'
#' @return A logic specifying whether \code{x} is of class
#'         \code{PredInactivationMCMC}
#'
#' @export
#'
is.PredInactivationMCMC <- function(x) inherits(x, "PredInactivationMCMC")

#------------------------------------------------------------------------------

## "summary" methods

#' Summary of a FitInactivation object
#'
#' @param object Instance of Fit Inactivation
#' @param ... ignored
#'
#' @export
#'
summary.FitInactivation <- function(object, ...) {

    summary(object$fit_results)

}

#' Summary of a FitInactivationMCMC object
#'
#' @param object Instance of FitInactivationMCMC
#' @param ... ignored
#'
#' @export
#'
summary.FitInactivationMCMC <- function(object, ...) {

    summary(object$modMCMC)

}

#' Summary of a IsoFitInactivation object
#'
#' @param object Instance of IsoFitInactivation
#' @param ... ignored
#'
#' @export
#'
summary.IsoFitInactivation <- function(object, ...) {

    summary(object$nls)

}




















