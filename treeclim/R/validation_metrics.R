## two commonly used metrics to validate reconstruction skills

##' @keywords internal
reduction_of_error <- function(observed, estimate, calibration) {
  1 - (sum((observed - estimate)^2)/
       (sum((observed - mean(calibration, na.rm = TRUE))^2)))
}

##' @keywords internal
coefficient_of_efficiency <- function(observed, estimate, verification) {
  1 - (sum((observed - estimate)^2)/
       (sum((observed - mean(verification, na.rm = TRUE))^2)))
}


