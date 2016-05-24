#' Function that downscales a coalescent model
#'
#' This function reduces the number of loci in all averaged loci by a
#' certain factor.
#' Non-averaged loci as created with \code{\link{locus_single}} are not
#' modified in any way. This function is primarily designed for jaatha,
#' and might be unsuitable for other purposes.
#'
#' @param model The model to downscale
#' @param scaling_factor The factor by which the number of loci are reduced.
#'   A value of 2 changes to numbers to half their value (rounded),
#'   a value of 3 to a third an so on.
#' @export
#' @examples
#' model <- coal_model(10, loci_number = 10) + locus_single(100)
#' model
#' # Group 1: 10 loci; group 2: 1 locus
#'
#' model <- scale_model(model, 3)
#' model
#' # Group 1: 3 loci; group 2: 1 locus
scale_model <- function(model, scaling_factor) {
  assert_that(is.model(model))
  assert_that(is.numeric(scaling_factor))
  assert_that(length(scaling_factor) == 1)
  model$scaling_factor <- scaling_factor
  model$id <- get_id()
  model
}
