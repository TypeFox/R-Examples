#' @importFrom R6 R6Class
mutation_class <- R6Class("mutation", inherit = feature_class,
  private = list(
    model = NA,
    fixed = FALSE,
    base_frequencies = NA,
    tstv_ratio = NA,
    gtr_rates = NA
  ),
  public = list(
    initialize = function(rate, model, base_frequencies,
                          tstv_ratio, gtr_rates, fixed) {
      private$rate <- private$add_parameter(rate, add_par = FALSE)

      assert_that(length(model) == 1)
      assert_that(any(model == c("IFS", "HKY", "GTR")))
      private$model <- model

      assert_that(length(fixed) == 1)
      assert_that(is.logical(fixed))
      private$fixed <- fixed

      if (model == "HKY") {
        if (is.na(tstv_ratio)) {
          stop("You need to specify tstv_ratio for the HKY mutation model")
        }
        assert_that(all(!is.na(base_frequencies)))
        assert_that(is.numeric(tstv_ratio))
        assert_that(length(tstv_ratio) == 1)
        private$tstv_ratio <- tstv_ratio

        if (any(is.na(base_frequencies))) {
          stop("missing base_frequencies for the HKY mutation model")
        }
        assert_that(all(!is.na(base_frequencies)))
        assert_that(is.numeric(base_frequencies))
        assert_that(length(base_frequencies) == 4)
        assert_that(sum(base_frequencies) == 1)
        private$base_frequencies <- base_frequencies
      }

      else if (model == "GTR") {
        if (any(is.na(gtr_rates))) {
          stop("You need to specify gtr_rates for the GTR mutation model")
        }
        assert_that(is.numeric(gtr_rates))
        assert_that(length(gtr_rates) == 6)
        private$gtr_rates <- gtr_rates
      }
    },
    get_model = function() private$model,
    get_base_frequencies = function() private$base_frequencies,
    get_tstv_ratio = function() private$tstv_ratio,
    get_gtr_rates = function() private$gtr_rates,
    get_fixed = function() private$fixed,
    print = function() {
      cat("Mutations with rate", print_par(paste0("par(", private$rate, ")")),
          "following a", private$model, "mutation model\n")
    }
  )
)

#' Feature: Mutation
#'
#' This feature adds mutations to a model. Mutations occur in the genomes
#' of the individuals with a given \code{rate}. The rate is per locus
#' for \link[=locus]{unlinked loci} and per trio for linked
#' \link[=locus_trio]{locus trios}. By default, the same mutation rate is used
#' for all loci, but it is possible to change this with \code{\link{par_variation}}
#' and \code{\link{par_zero_inflation}}.
#'
#' @param rate The mutation rate. Can be a numeric or a \code{\link{parameter}}.
#'        The rate is specified as \eqn{4 * N0 * mu}, where \eqn{mu} is the
#'        mutation rate per locus.
#' @param fixed_number If set to \code{TRUE}, the number of mutations on each
#'   locus will always be exactly equal to the rate, rather than happening with
#'   a rate along the ancestral tree.
#' @param model The mutation model you want to use.
#'   Can be either 'IFS' (default), 'HKY' or 'GTR'. Refer to the mutation model
#'   section for detailed information.
#' @param tstv_ratio The ratio of transitions to transversions used in the 'HKY'
#'   muation model.
#' @param base_frequencies The equilibrium frequencies of the four bases used in
#'   the 'HKY' mutation model. Must be a numeric vector of length four, with the
#'   values for A, C, G and T, in that order.
#' @param gtr_rates The rates for the six amino acid substitutions used in the
#'   'GTR' model. Must be a numeric vector of length six.
#'   Order: A<->C, A<->G, A<->T, C<->G, C<->T, G<->T.
#' @return The feature, which can be added to a model using `+`.
#' @export
#' @seealso For using rates that variate between the loci in a model:
#'   \code{\link{par_variation}}, \code{\link{par_zero_inflation}}
#' @seealso For adding recombination: \code{\link{feat_recombination}}.
#' @family features
#'
#' @section Mutation Models:
#' The infinite sites mutation (\strong{IFS}) model is a frequently used simplification
#' in population genetics. It assumes that each locus consists of infinitely
#' many sites at which mutations can occur, and each mutation hits a new site.
#' Consequently, there are no back-mutations with this model. It does not
#' generate DNA sequences, but rather only 0/1 coded data, were 0 denotes the
#' ancestral state of the site, and 1 the derived state created by a mutation.
#'
#' The other mutation models are finite site models that generate more realistic
#' sequences.
#'
#' The Hasegawa, Kishino and Yano (\strong{HKY}) model (Hasegawa et al., 1985) allows
#' for a different rate of transitions and transversions (tstv_ratio)
#' and unequal
#' frequencies of the four nucleotides (base_frequencies).
#'
#' The general reversible process (\strong{GTR}) model (e.g. Yang, 1994) is more general
#' than the HKY model and allows to define the rates for each
#' type of substitution. The rates are assumed to be symmetric
#' (e.g., the rate for T to G is equal to the one for G to T).
#'
#' @examples
#' # A model with a constant mutation rate of 5:
#' model <- coal_model(5, 1) + feat_mutation(5) + sumstat_seg_sites()
#' simulate(model)
#'
#' # A model with 7 mutations per locus:
#' model <- coal_model(5, 1) + feat_mutation(7, fixed = TRUE) + sumstat_seg_sites()
#' simulate(model)
#'
#' # A model using the HKY model:
#' model <- coal_model(c(10, 1), 2) +
#'  feat_mutation(7.5, model = "HKY", tstv_ratio = 2,
#'                base_frequencies = c(.25, .25, .25, .25)) +
#'  feat_outgroup(2) +
#'  feat_pop_merge(1.0, 2, 1) +
#'  sumstat_seg_sites()
#'  \dontrun{simulate(model)}
#'
#' # A model using the GTR model:
#' model <- coal_model(c(10, 1), 1, 25) +
#'  feat_mutation(7.5, model = "GTR",
#'                gtr_rates = c(1, 1, 1, 1, 1, 1) / 6) +
#'  feat_outgroup(2) +
#'  feat_pop_merge(1.0, 2, 1) +
#'  sumstat_dna()
#'  \dontrun{simulate(model)$dna}
feat_mutation <- function(rate,
                          model = "IFS",
                          base_frequencies = NA,
                          tstv_ratio = NA,
                          gtr_rates = NA,
                          fixed_number = FALSE) {

  mutation_class$new(rate, model, base_frequencies,
                     tstv_ratio, gtr_rates, fixed_number)
}

is_feat_mutation <- function(feat) any("mutation" == class(feat))

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.mutation <- function(feature, model) {
  if (feature$get_model() != "IFS") stop("Unsupported mutation model")
  if (feature$get_fixed()) {
    paste0("-s', par(", feature$get_rate(), "), '")
  } else {
    paste0("-t', par(", feature$get_rate(), "), '")
  }
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.mutation <- conv_to_ms_arg.mutation

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.mutation <- function(feature, model) {
  if (feature$get_fixed()) {
    stop("scrm does not support simulating a fixed number of mutations",
         call. = FALSE)
  }
  conv_to_ms_arg.mutation(feature, model)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.mutation <- function(feature, model) {
  if (feature$get_model() == "IFS") {
    stop("seq-gen can not simulate an IFS model", call. = FALSE)
  }
  if (feature$get_fixed()) {
    stop("seq-gen can not simulate a fixed number of mutations", call. = FALSE)
  }
  if (feature$get_model() == "GTR") {
    rates <- paste("-r", paste(feature$get_gtr_rates(), collapse = " "))
  } else {
    rates <- paste("-f", paste(feature$get_base_frequencies(), collapse = " "),
                   "-t", feature$get_tstv_ratio())
  }
  paste0("-m", feature$get_model(), " ",
         rates, " ",
         "-l', locus_length, '",
         "-s', par(", feature$get_rate(), " / locus_length), '",
         "-p', locus_length + 1, '",
         "-z', par(seed), '-q")
}
