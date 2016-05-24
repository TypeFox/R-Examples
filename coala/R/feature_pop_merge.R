pop_merge_class <- R6Class("pop_merge", inherit = feature_class,
  public = list(
    initialize = function(time, pop_from, pop_to) {
      private$time <- private$add_parameter(time)
      private$set_population(c(from = pop_from, to = pop_to), 2)
    },
    print = function() {
      cat("Merge of pop", private$population[1],
          "into pop", private$population[2],
          "at time", print_par(self$get_time()), "\n")
    }
  )
)


#' Feature: Population Merge
#'
#' Backwards in time, this feature merges one population into another.
#' Forwards in time, this corresponds to a speciation event.
#'
#' In addition to the merge, the growth rate of and all migration rates from
#' the source population will be set to 0 at the time of the merge to mimic
#' a speciation event forwards in time. Technically, \code{pop_source} is
#' still present in the model, but not used unless migration to the population
#' is manually enabled.
#'
#' @param pop_source The population from which all lines are moved.
#'        This is the newly created population in the speciation event.
#' @param pop_target The population to which the lines are moved. This is
#'        the population in which the speciation event occurs.
#' @param time The time at which the merge occurs.
#' @export
#' @family features
#' @examples
#' # Two population which merge after 0.5 time units:
#' model <- coal_model(c(25,25), 1) +
#'   feat_pop_merge(0.5, 2, 1) +
#'   feat_mutation(5) +
#'   sumstat_tajimas_d()
#' simulate(model)
#'
#' # An standard isolation-with-migration model:
#' model_iwm <- model +
#'   feat_migration(.75, symmetric = TRUE)
#' simulate(model_iwm)
feat_pop_merge <- function(time, pop_source, pop_target) {
  pop_merge_class$new(time, pop_source, pop_target)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.pop_merge <- function(feature, model) {
  paste0("-ej', ", feature$get_time(), ", ",
         feature$get_population()[1], ", ",
         feature$get_population()[2], ", '")
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.pop_merge <- conv_to_ms_arg.pop_merge

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.pop_merge <- conv_to_ms_arg.pop_merge

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.pop_merge <- ignore_par
