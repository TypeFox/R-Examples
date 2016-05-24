growth_class <- R6Class("growth", inherit = feature_class,
  public = list(
    print = function() {
      cat("Exponential growth/decline with rate", print_par(private$rate),
          "in population", self$get_population(),
          "starting at time", print_par(self$get_time()), "\n")
    }
  )
)

#' Feature: Exponential population size growth/decline
#'
#' This feature changes the growth factor of a population at given
#' point in time. This factor applies to the time interval further
#' into the past from this point.
#'
#' The population size changes by a factor \eqn{exp(-\alpha*t)}, where
#' \eqn{\alpha} is the growth parameter and \eqn{t} is the time since
#' the growth has started. For positive alpha, the population will decline
#' backwards in time or grow forwards in time. For a negative value of
#' \eqn{\alpha} it will decline (forward in time).
#'
#' @param rate The growth rate. Can be a numeric or a \code{\link{parameter}}.
#'        See \code{Details} for an explanation how the rate affects the
#'        population size.
#' @param population The population which growths/declines. Can be
#'          "all" for all populations, or the number of one population.
#' @param time The time at which the growth rate is changed. Can also be
#'        a \code{\link{parameter}}.
#' @export
#' @seealso For instantaneous population size
#'          changes: \code{\link{feat_size_change}}
#' @family features
#' @examples
#' # Simulate a haploid population that has been expanding for
#' # the last 2*Ne generations
#' model <- coal_model(10, 1) +
#'   feat_growth(5, time = 0) +
#'   feat_growth(0, time = 1) +
#'   feat_mutation(10) +
#'   sumstat_sfs()
#' simulate(model)
feat_growth <- function(rate, population = "all", time="0") {
  growth_class$new(rate, population, time)
}


#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.growth <- function(feature, model) {
  all_pops <- feature$get_population() == "all" ||
    (feature$get_population() == 1 && length(get_populations(model)) == 1)
  present <- feature$get_time() == "par(0)"

  if (present) {
    if (all_pops) cmd <- "-G"
    else cmd <- "-g"
  } else {
    if (all_pops) cmd <- "-eG"
    else cmd <- "-eg"
  }

  paste0(cmd, "', ",
         ifelse(present, "", paste(feature$get_time(), ", ")),
         ifelse(all_pops, "", paste(feature$get_population(), ", ")),
         feature$get_rate(), ", '")
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.growth <- conv_to_ms_arg.growth

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.growth <- conv_to_ms_arg.growth

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.growth <- ignore_par
