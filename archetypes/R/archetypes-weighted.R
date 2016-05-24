#' @include archetypes-class.R
#' @include archetypes-kit.R
#' @include archetypes-kit-blocks.R
{}



#' Weighted archetypes
#'
#' @inheritParams archetypes
#' @param weights Data weights matrix.
#' @param familyBlocks Exchange predefined family blocks; see
#'   \code{\link{archetypesFamily}}.
#'
#' @return An object of class \code{weightedArchetypes} and
#'   \code{\link{as.archetypes}}.
#'
#' @family archetypes
#'
#' @export
weightedArchetypes <- function(data, k, weights = NULL,
                               familyBlocks = list(), ...) {

  family <- do.call(archetypesFamily, c(list('weighted'), familyBlocks))

  archetypes(data, k, weights = weights, family = family, ...)
}



.weighted.archetypesFamily <- function() {
  f <- .original.archetypesFamily()
  f$class <- 'weightedArchetypes'
  f$globweightfn <- center.globweightfn
  f
}



setOldClass(c("weightedArchetypes", "archetypes"))
