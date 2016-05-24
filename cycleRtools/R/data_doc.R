#' Example cycling data.
#'
#' @description Formatted cycling data from a Garmin head-unit. Imported via
#' \code{read_fit("file_path", format = TRUE, CP = 310, sRPE = 7)}.
#'
#' \code{"ridedata"} is a typical group ride. \code{"intervaldata"} is a session
#' (of sorts) that included two efforts and a cafe stop. The latter is included
#' to demonstrate the use of \code{\link{interval_detect}}.
#'
#' @format An object of class \code{c("cycleRdata", "data.frame")}, and
#'   additional attributes of \code{CP = 300} & \code{sRPE = 7}. The latter are
#'   used by several methods in this package. See \code{\link{cycleRdata}} for a
#'   description of columns.
#'
#' @seealso \code{\link{cycleRdata}}.
#'
#' @name ride_examples
NULL

#' @rdname ride_examples
"ridedata"
#' @rdname ride_examples
"intervaldata"

#' Power-time profile.
#'
#' An example power profile; i.e. best mean powers for periods of 30 seconds
#' through to 1 hour, in increments of 10 seconds.
#'
#' @format a \code{data.frame} with two columns: time (seconds) and power
#'   (Watts), respectively.
"Pt_prof"
