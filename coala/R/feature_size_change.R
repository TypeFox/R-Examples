size_change_class <- R6Class("size_change", inherit = feature_class,
  public = list(
    print = function() {
      cat("Instanious size change to", print_par(private$rate),
          "* N0 in population", self$get_population(),
          "at time", print_par(self$get_time()), "\n")
    }
  )
)


#' Feature: Instantaneous Size Change
#'
#' This feature changes the effective population size of one
#' population. The change is performed at a given time point
#' and applies to the time interval further on into
#' the past from this point. The population size is set to a
#' fraction of N0.
#'
#' @param new_size A \code{\link{parameter}} giving the new size of the
#'   population, as a fraction of N0.
#' @param population The number of the population whichs size changes.
#'          Can also be set to "all". Then the size changes applies to all
#'          populations.
#' @param time The time at which the population's size is changed.
#' @return The feature, which can be added to a model using `+`.
#' @export
#' @seealso For continuous size changes over time: \code{\link{feat_growth}}.
#' @family features
#' @examples
#' # A model with one smaller population:
#' model <- coal_model(c(20, 5), 3) +
#'   feat_size_change(.1, population = 2) +
#'   feat_mutation(1.0) +
#'   feat_migration(0.5, 2, 1) +
#'   sumstat_sfs()
#' simulate(model)
#'
#' # A model of one population that experienced a bottleneck:
#' model <- coal_model(10, 1) +
#'   feat_size_change(0.1, time = 0.3) +
#'   feat_size_change(1.0, time = 0.5) +
#'   feat_mutation(20) +
#'   sumstat_sfs()
#' simulate(model)
feat_size_change <- function(new_size, population = 1, time = "0") {
  size_change_class$new(new_size, population, time)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.size_change <- function(feature, model) {
  all_pops <- feature$get_population() == "all" ||
    (feature$get_population() == 1 && length(get_populations(model)) == 1)
  present <- feature$get_time() == "par(0)"

  if (present) {
    if (all_pops) cmd <- "-N"
    else cmd <- "-n"
  } else {
    if (all_pops) cmd <- "-eN"
    else cmd <- "-en"
  }

  paste0(cmd, "', ",
         ifelse(present, "", paste(feature$get_time(), ", ")),
         ifelse(all_pops, "", paste(feature$get_population(), ", ")),
         feature$get_rate(), ", '")
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.size_change <- conv_to_ms_arg.size_change

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.size_change <- conv_to_ms_arg.size_change

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.size_change <- ignore_par
