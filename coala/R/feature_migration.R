migration_class <- R6Class("migration", inherit = feature_class,
  private = list(rate = NA),
  public = list(
    initialize = function(rate, pop_from, pop_to, time, symmetric = FALSE) {
      super$initialize(rate = rate, time = time)

      if (symmetric) {
        private$population <- "all"
      } else {
        private$set_population(c(from = pop_from, to = pop_to), 2)
      }
    },
    print = function() {
      if (all(private$population == "all")) {
        cat("Symmetric migration")
      } else {
        cat("Migration from pop", private$population,
            "to pop", private$pop_to)
      }
      cat(" with rate", print_par(private$rate),
          "starting at time", print_par(self$get_time()), "\n")
    }
  )
)


#' Feature: Migration/Gene Flow
#'
#' This feature changes the migration rates at a given time point.
#' Per default, no migration between the population occurs, which corresponds
#' to a \code{rate} of \code{0}. Set it to a value greater than zero to
#' enable migration from one population to another.
#'
#' When looking forward in time, a fraction of \code{pop_to} that is replaced
#' by migrants from \code{pop_from} each generation (see \code{rate}). When
#' looking backwards in time, ancestral lines in \code{pop_to} move to
#' \code{pop_from} with the given rate.
#'
#' @param rate The migration rate. Can be a numeric or a
#'        \code{\link{parameter}}. The rate is specified as
#'        \eqn{4 * N0 * m}, where \eqn{m} is the fraction of
#'        \code{pop_to} that is replaced by migrants
#'        from \code{pop_from} each generation (in forward time).
#' @param pop_from The population from which the individuals leave.
#' @param pop_to The population to which the individuals move.
#' @param symmetric Use the rate for all pairs of populations.
#' @param time The time point at which the migration with the migration
#'        rate is set. The rate applies to the time past warts
#'        of the time point, until it is changed again.
#' @export
#' @family features
#'
#' @examples
#' # Asymmetric migration between two populations:
#' model <- coal_model(c(5, 5), 10) +
#'   feat_migration(0.5, 1, 2) +
#'   feat_migration(1.0, 2, 1) +
#'   feat_mutation(5) +
#'   sumstat_sfs()
#' simulate(model)
#'
#' # Three populations that exchange migrations with equal
#' # rates at times more than 0.5 time units in the past:
#' model <- coal_model(c(3, 4, 5), 2) +
#'   feat_migration(1.2, symmetric = TRUE, time = 0.5) +
#'   feat_mutation(5) +
#'   sumstat_sfs()
#' simulate(model)
feat_migration <- function(rate, pop_from = NULL, pop_to = NULL,
                           symmetric = FALSE, time = "0") {
  if (symmetric) {
    return(migration_class$new(rate, time = time, symmetric = TRUE))
  } else {
    return(migration_class$new(rate, pop_from, pop_to, time))
  }
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.migration <- function(feature, model) {
  if (all(feature$get_population() == "all")) {
    return( paste0("-eM', ", feature$get_time(), ", ",
                    feature$get_rate(), ", '"))
  }
  paste0("-em', ", feature$get_time(), ", ",
         feature$get_population()[2], ", ",
         feature$get_population()[1], ", ",
         feature$get_rate(), ", '")
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.migration <- conv_to_ms_arg.migration

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.migration <- conv_to_ms_arg.migration

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.migration <- ignore_par
