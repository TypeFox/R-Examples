unphased_class <- R6Class("unphased", inherit = feature_class,
  private = list(ploidy = NA, samples_per_ind = NA),
  public = list(
    initialize = function(samples_per_ind) {
      assert_that(is.number(samples_per_ind))
      private$samples_per_ind <- samples_per_ind
    },
    check = function(model) {
      if (self$get_samples_per_ind() > get_ploidy(model)) {
        stop("samples_per_ind needs to be lower or equal to the ploidy",
             call. = FALSE)
      }
      invisible(NULL)
    },
    get_samples_per_ind = function() private$samples_per_ind,
    print = function() {
      cat("Unphasing of the simulated chromosomes into",
          private$samples_per_ind, "pseudo-chromosomes per individual\n")
    }
  )
)

#' Feature: Unphased Sequences
#'
#' This simulates unphased data by randomly mixing the sites within one
#' individual. Each position is randomly taken from a phased chromosome..
#'
#' For each individual, \code{ploidy} chromosomes are simulated, and
#' \code{samples_per_ind} pseudo-chromosomes are created of these.
#'
#' @param samples_per_ind The number of pseudo-chromosomes that are created
#'   from the phased chromosomes for each individual.
#' @return The feature, which can be added to a model using `+`.
#' @export
#' @family features
#' @examples
#' # Simulate unphased data in a diploid population
#' model <- coal_model(10, 1, ploidy = 2) +
#'   feat_mutation(10) +
#'   feat_unphased(2) +
#'   sumstat_seg_sites()
#' simulate(model)
#'
#' # The same as before, but return only one chromosome for
#' # each individual:
#' model <- coal_model(10, 1, ploidy = 2) +
#'   feat_mutation(10) +
#'   feat_unphased(1) +
#'   sumstat_seg_sites()
#' simulate(model)
feat_unphased <- function(samples_per_ind) {
  unphased_class$new(samples_per_ind)
}


is_feat_unphased <- function(feat) any("unphased" == class(feat))


get_feature_unphased <- function(model) {
  mask <- vapply(model$features, is_feat_unphased, logical(1))
  if (!any(mask)) return(NULL)
  if (sum(mask) > 1) stop("multiple unphasings are not supported")
  model$features[mask][[1]]
}


get_samples_per_ind <- function(model) {
  feat <- get_feature_unphased(model)
  if (is.null(feat)) return(get_ploidy(model))
  as.integer(feat$get_samples_per_ind())
}


is_unphased <- function(model) !is.null(get_feature_unphased(model))

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.unphased <- ignore_par

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.unphased <- ignore_par

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.unphased <- ignore_par

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.unphased <- ignore_par
