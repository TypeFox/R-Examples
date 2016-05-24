selection_class <- R6Class("selection", inherit = feature_class,
  private = list(
    strength_AA = 0,
    strength_Aa = 0,
    strength_aa = 0,
    additive = FALSE,
    force_keep = TRUE,
    start_frequency = 0.05,
    Ne = 10000,
    position = 0.5,
    start = TRUE
  ),
  public = list(
    initialize = function(strength_AA, strength_Aa, strength_aa,
                          population, time, additive = FALSE,
                          start, start_frequency, Ne, position, force_keep) {

      assert_that(is.logical(additive) && length(additive) == 1)
      private$additive <- additive

      assert_that(is.logical(start) && length(start) == 1)
      private$start <- start

      assert_that(is.numeric(start_frequency))
      private$start_frequency <- start_frequency

      assert_that(is.number(Ne))
      private$Ne <- Ne

      assert_that(is.numeric(position) && length(position) == 1)
      private$position <- position

      assert_that(is.logical(force_keep) && length(force_keep) == 1)
      private$force_keep <- force_keep

      private$strength_Aa <- private$add_parameter(strength_Aa)
      if (!additive) {
        private$strength_AA <- private$add_parameter(strength_AA)
        private$strength_aa <- private$add_parameter(strength_aa)
      }

      private$time <- private$add_parameter(time)
      private$set_population(population)
    },
    get_strength_AA = function() private$strength_AA,
    get_strength_Aa = function() private$strength_Aa,
    get_strength_aa = function() private$strength_aa,
    is_additive = function() private$additive,
    get_start = function() private$start,
    get_force_keep = function() private$force_keep,
    get_Ne = function() private$Ne,
    get_start_freq = function() private$start_frequency,
    get_position = function() private$position,

    print = function() {
      if (private$additive) {
        cat("Additive selection with strength ", print_par(private$strength_AA))
      } else {
        cat("Selection with strength ", print_par(private$strength_AA),
            "(AA) and ", print_par(private$strength_Aa), "(Aa)")
      }
      cat(" in population", self$get_population(),
          "starting at time", print_par(self$get_time()), "\n")
    }
  )
)

#' Feature: Selection
#'
#' This feature adds selection to a model. Only one site per locus can be under
#' selection. Using this feature requires that \code{msms} is installed, see
#' \code{\link{activate_msms}}.
#'
#' @param population The population in which the allele is selected. Can either
#'   be \code{all} for all population, or the number of a population.
#' @param time The time at which the selection starts if \code{start == TRUE}
#'   (looking forwards in time), or the time at which the selection strength
#'   changes if \code{start == FALSE}. The new strength applies for the time
#'   period further into the past in this case.
#' @param strength_AA The selection strength for the selected homozygote.
#'   The parameter is valid for the chosen population and the time further
#'   past-wards from either time 0 on if \code{start = TRUE}, or from \code{time}
#'   onwards. The same applies for \code{strength_Aa}, \code{strength_aa} and
#'   \code{strength_A}.
#' @param strength_Aa The selection strength for the heterozygote.
#' @param strength_aa The selection strength for the recessive homozygote.
#' @param strength_A This sets the strength for the selected allele in a
#'   haploid model or a diploid model with additive selection.
#'   \code{strength_AA}, \code{strength_Aa}, \code{strength_aa}
#'   are ignored when this is argument is given.
#' @param start Whether selection should start at this time point. At this time,
#'   the selected allele is introduced in the population with an initial
#'   starting frequency. This must be set to \code{TRUE} for exactly one
#'   selection feature in the model. The values of \code{start_frequency},
#'   \code{Ne}, \code{position} and \code{force_keep} are used for the model.
#'   You can add additional selection feature to the model to set the
#'   selection strength for more demes or change it at different time points,
#'   but these need to have \code{start = FALSE}.
#' @param start_frequency The start frequency at which the selected allele is
#'   introduced at \code{time}. If the model has multiple population, this
#'   can either be a numeric vector that contains the initial frequency for each
#'   population or a single number. In the latter case, the value is used for
#'   all population specified with \code{populations}, and 0 is used for all
#'   other populations.
#' @param Ne The effective population size that is used for the forward
#'   simulations.
#' @param position The position of the selected site, relative to the
#'   simulated sequence. Values between 0 and 1 are within the simulated area,
#'   while smaller values are to the left of it and larger ones to the right.
#' @param force_keep Whether to restart simulations in which the selected goes to
#'   extinction or not.
#'
#' @export
#' @seealso For using rates that variate between the loci in a model:
#'   \code{\link{par_variation}}, \code{\link{par_zero_inflation}}
#' @seealso For summary statistics that are sensitive for selection:
#'   \code{\link{sumstat_tajimas_d}}, \code{\link{sumstat_ihh}},
#'   \code{\link{sumstat_omega}}, \code{\link{sumstat_mcmf}}
#' @family features
#' @examples
#' # Positive additive selection in population 2:
#' model <- coal_model(c(10, 13), 1, 10000) +
#'   feat_pop_merge(.5, 2, 1) +
#'   feat_selection(strength_A = 1000,
#'                  population = 2,
#'                  time = par_named("tau")) +
#'   feat_mutation(100) +
#'   feat_recombination(10) +
#'   sumstat_tajimas_d(population = 2)
#' \dontrun{simulate(model, pars = c(tau = 0.03))}
feat_selection <- function(strength_AA = 0,
                           strength_Aa = 0,
                           strength_aa = 0,
                           strength_A = NULL,
                           population = "all",
                           time,
                           start = TRUE,
                           start_frequency = 0.0005,
                           Ne = 10000,
                           position = 0.5,
                           force_keep = TRUE) {

  if (!is.null(strength_A)) {
    return(selection_class$new(strength_Aa = strength_A,
                               population = population,
                               time = time,
                               additive = TRUE,
                               start = start,
                               start_frequency = start_frequency,
                               Ne = Ne,
                               position = position,
                               force_keep = force_keep))
  }

  selection_class$new(strength_AA, strength_Aa, strength_aa, population, time,
                      FALSE, start, start_frequency, Ne, position, force_keep)
}


#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.selection <- function(feature, model) {
  stop("selection is not supported", call. = FALSE)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.selection <- conv_to_ms_arg.selection

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.selection <- function(feature, model) {

  # Determine general options if we are starting selection
  if (feature$get_start()) {
    n_pop <- length(get_populations(model))

    # Determine initial frequency of selected allele
    if (length(feature$get_start_freq()) == 1) {
      if (feature$get_population() == "all") {
        start_freq <- rep(feature$get_start_freq(), n_pop)
      } else {
        start_freq <- rep(0, n_pop)
        start_freq[feature$get_population()] <- feature$get_start_freq()
      }
    } else {
      start_freq <- feature$get_start_freq()
    }

    start_cmd <- paste0("-SI', ", feature$get_time(), ", '",
                        n_pop, " " , paste(start_freq, collapse = " "), " ",
                        "-Sp', ", feature$get_position(), ", '",
                        "-N', ", feature$get_Ne(), ", '",
                        ifelse(feature$get_force_keep(), "-SForceKeep ", ""))
  } else {
    start_cmd <- ""
  }

  if (feature$is_additive()) {
    if (feature$get_population() == "all") {
      strength <- paste0("-SA',", feature$get_strength_Aa(), ", '")
    } else {
      strength <- paste0("-Sc',",
                         ifelse(feature$get_start(), 0, feature$get_time()), ", ", #nolint
                         feature$get_population(), ", ",
                         feature$get_strength_Aa(), ", '")
    }
  } else {
    if (feature$get_population() == "all") {
      strength <- paste0("-SAA',", feature$get_strength_AA(), ", '",
                         "-SAa',", feature$get_strength_Aa(), ", '",
                         "-Saa',", feature$get_strength_aa(), ", '")
    } else {
      strength <- paste0("-Sc',",
                         ifelse(feature$get_start(), 0, feature$get_time()), ", ", #nolint
                         feature$get_population(), ", ",
                         feature$get_strength_AA(), ", ",
                         feature$get_strength_Aa(), ", ",
                         feature$get_strength_aa(), ", '")
    }
  }
  paste0(start_cmd, strength)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.selection <- ignore_par
