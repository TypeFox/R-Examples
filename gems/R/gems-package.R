#' gems: Generalized Multistate Simulation Model
#'
#' Simulate and analyze multistate models with 
#' general hazard functions. gems provides functionality for the preparation 
#' of hazard functions and parameters, simulation from a general multistate 
#' model and predicting future events. The multistate model is not required 
#' to be a Markov model and may take the history of previous events into 
#' account. In the basic version, it allows to simulate from 
#' transition-specific hazard function, whose parameters are multivariable 
#' normally distributed.
#'
#' @docType package
#' @name gems
#' @references 
#' Nello Blaser, Luisa Salazar Vizcaya, Janne Estill, Cindy Zahnd, 
#' Bindu Kalesan, Matthias Egger, Olivia Keiser, Thomas Gsponer (2015). gems: 
#' An R Package for Simulating from Disease Progression Models. Journal of 
#' Statistical Software, 64(10), 1-22. URL http://www.jstatsoft.org/v64/i10/.
NULL

#' tavi data set
#' 
#' The simulated data set for each patient contains data for kidney injuries,
#' bleeding complications and the combined endpoint of stroke or death. The
#' data was simulated from the original data following the steps described in
#' the package vignette.
#' 
#' 
#' @name tavi
#' @docType data
#' @format A data frame with 194 observations on the following 7 variables.
#' \describe{
#' \item{id}{a character vector that contains the patient
#' id's} 
#' \item{kidney}{a numeric vector; indicator variable that show
#' if an event has occurred} 
#' \item{kidney.dur}{a numeric vector; times
#' at which the events occurred or the patients were censored}
#' \item{bleeding}{a numeric vector; indicator variable that show if an
#' event has occurred} 
#' \item{bleeding.dur}{a numeric vector; times at
#' which the events occurred or the patients were censored}
#' \item{death}{a numeric vector; indicator variable that show if an
#' event has occurred} 
#' \item{death.dur}{a numeric vector; times at
#' which the events occurred or the patients were censored} 
#' }
#' @keywords datasets
#' @examples
#' 
#' head(data(tavi))
#' 
#' @references 
#' Nello Blaser, Luisa Salazar Vizcaya, Janne Estill, Cindy Zahnd, 
#' Bindu Kalesan, Matthias Egger, Olivia Keiser, Thomas Gsponer (2015). gems: 
#' An R Package for Simulating from Disease Progression Models. Journal of 
#' Statistical Software, 64(10), 1-22. URL http://www.jstatsoft.org/v64/i10/.
NULL