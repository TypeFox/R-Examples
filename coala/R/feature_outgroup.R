outgroup_class <- R6Class("outgroup", inherit = feature_class,
  public = list(
    initialize = function(population) {
      private$set_population(population)
    },
    print = function() {
      cat("Outgroup: Population", private$population, "\n")
    }
  )
)

#' Feature: Outgroup
#'
#' This feature declares an existing population as outgroup. Outgroups are used
#' to determine the ancestral allele in finite sites simulations and are required
#' there. All individuals of the outgroup are ignored when calculating summary
#' statistics. If the outgroup consists of multiple individuals, all positions
#' where the individuals have different alleles are ignored.
#'
#' @param population The population that is marked as outgroup.
#' @export
#' @family features
#' @examples
#' # A simple finite sites model
#' model <- coal_model(c(4, 6, 1), 2, 10) +
#'    feat_outgroup(3) +
#'    feat_pop_merge(0.5, 2, 1) +
#'    feat_pop_merge(2, 3, 1) +
#'    feat_mutation(5, model = "GTR", gtr_rates = 1:6)
feat_outgroup <- function(population) {
  outgroup_class$new(population)
}

is_feat_outgroup <- function(feat) inherits(feat, "outgroup")


#' @describeIn get_features Returns the population that is marked as outgroup
#' @export
get_outgroup <- function(model) {
  mask <- vapply(model$features, is_feat_outgroup, logical(1))
  if (sum(mask) != 1) return(NA)
  model$features[mask][[1]]$get_population()
}


#' @describeIn get_features Returns the number of samples in the outgroup
#' @export
get_outgroup_size <- function(model, for_sim = FALSE) {
  outgroup <- get_outgroup(model)
  if (is.na(outgroup)) return(0)
  outgroup_size <- get_sample_size(model, for_sim)[outgroup]
  outgroup_size
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.outgroup <- function(feature, model) {
  stop("ms does not support outgroups", call. = FALSE)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.outgroup <- function(feature, model) {
  stop("msms does not support outgroups", call. = FALSE)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.outgroup <- function(feature, model) {
  stop("scrm does not support outgroups", call. = FALSE)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.outgroup <- ignore_par
