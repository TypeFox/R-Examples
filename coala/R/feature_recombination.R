recombination_class <- R6Class("recombination", inherit = feature_class,
  public = list(
    print = function() {
      cat("Recombination with rate", print_par(private$rate), "\n")
    }
  )
)

#' Feature: Recombination
#'
#' This feature adds intra-locus recombination to a model.  The rate is per locus
#' for \link[=locus]{unlinked loci} and per trio for linked
#' \link[=locus_trio]{locus trios}. By default, the same recombination rate is used
#' for all loci, but it is possible to change this with \code{\link{par_variation}}
#' and \code{\link{par_zero_inflation}}. Coala assumes that recombination events
#' can occur between all neighboring bases.
#'
#' @param rate The recombination rate. Can be a numeric or a
#'        \code{\link{parameter}}. The rate is equal to
#'        \eqn{4*N0*r}, where \eqn{r} is the probability that a
#'        recombination event within the locus occurs in one generation.
#' @return The feature, which can be added to a model using `+`.
#' @export
#' @family features
#' @seealso For adding recombination: \code{\link{feat_mutation}}.
#'
#' @examples
#' # Simulate a 1.5kb sequence for 10 individuals with recombination:
#' model <- coal_model(10, 2, 1500) +
#'   feat_recombination(1.5) +
#'   feat_mutation(5) +
#'   sumstat_sfs()
#' simulate(model)
feat_recombination <- function(rate) {
  recombination_class$new(rate)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.recombination <- function(feature, model) {
  paste0("-r', ", feature$get_rate(), ", par(locus_length), '")
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.recombination <- conv_to_ms_arg.recombination

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.recombination <- conv_to_ms_arg.recombination

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.recombination <- ignore_par
