#' @include archetypes-class.R
#' @include archetypes-kit.R
#' @include archetypes-kit-blocks.R
{}



#' Robust archetypes
#'
#' @inheritParams archetypes
#' @param familyBlocks Exchange predefined family blocks; see
#'   \code{\link{archetypesFamily}}.
#'
#' @return An object of class \code{robustArchetypes} and
#'   \code{\link{as.archetypes}}.
#'
#' @family archetypes
#'
#' @export
robustArchetypes <- function(data, k, familyBlocks = list(), ...) {

  family <- do.call(archetypesFamily, c(list('robust'), familyBlocks))

  archetypes(data, k, family = family, ...)
}



.robust.archetypesFamily <- function() {
  f <- .original.archetypesFamily()
  f$class <- 'robustArchetypes'
  f$weightfn <- center.weightfn
  f$reweightsfn <- bisquare0.reweightsfn
  f
}



setOldClass(c("robustArchetypes", "archetypes"))
