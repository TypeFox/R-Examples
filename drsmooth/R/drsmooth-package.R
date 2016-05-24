#' drsmooth-package:  Dose-response Modeling with Smoothing Splines
#' 
#' drsmooth provides tools for assessing the shape of a dose-response curve
#' by testing linearity and non-linearity at user-defined cut-offs.
#' It also provides two methods of estimating a threshold dose,
#' or the dose at which the dose-response function transitions to significantly
#' increasing:  bi-linear (based on pkg:segmented)
#' and smoothed with splines (based on pkg:mgcv).
#' 
#' v.1.9.0 introduces spline-based dose-response modeling on dichotomous
#' data, with examples using the included DIdata.  See NEWS for details,
#' as well as the help files for the exported functions itemized below.
#' 
#' @docType package
#' @name drsmooth-package
#' @seealso drsmooth
#'
#' @section drsmooth functions:
#' There are 8 user-initiated functions in this package; see the help pages
#' for documentation of each.
#'
#' prelimstats() executes up to 5 tests of homogeneity & normality.
#'
#' noel() executes up to 5 tests to determine the no-observed-effect-level.
#'
#' nlaad(), lbcd(), nlbcd() test linearity across all doses, linearity below a
#' defined cut-off dose, and non-linearity below a defined cut-off dose,
#' respectively.
#' 
#' spline.plot() only prints the smoothed dose-response curve.
#'
#' segment() returns a two-segment linear dose-response model, by imposing a zero
#'  slope for the initial (left) segment and detecting one breakpoint where the
#' dose-response relation changes to a positive slope (if such a breakpoint exists).
#'
#' drsmooth() generates a spline model with the input dose and response,
#' plots the spline-estimated dose-response function with its upper and lower
#' 95 percent confidence bounds along with the actual data,
#' and returns key metrics related to the dose-response function.
#' 
#' expand () expands summarized dichotomous data into the format expected
#' by drsmooth(), if needed.

NULL