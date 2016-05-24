#' Returns a character vector of possible function tags.
#'
#' Test function are frequently distinguished by characteristic high-level properties,
#' e.g., unimodal or multimodal, continuous or discontinuous, separable or non-separable.
#' The \pkg{smoof} package offers the possibility to associate a set of properties,
#' termed \dQuote{tags} to a \code{smoof_function}. This helper function returns
#' a character vector of all possible tags.
#'
#' @return [\code{character}]
#' @export
getAvailableTags = function() {
  c("unimodal", "multimodal",
    "separable", "non-separable",
    "convex", "non-convex",
    "continuous", "discontinuous",
    "scalable", "non-scalable",
    "differentiable", "non-differentiable",
    "low-conditioned", "moderate-conditioned", "highly-conditioned",
    "adequate-global-structure", "weak-global-structure",
    "single-objective", "multi-objective"
  )
}
