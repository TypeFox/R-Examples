#' @importFrom R6 R6Class
stat_pi_class <- R6Class("stat_pi", inherit = sumstat_class,
  private = list(
    population = NULL,
    req_segsites = TRUE
  ),
  public = list(
    initialize = function(name, population, transformation) {
      assert_that(length(population) == 1)
      private$population <- population
      super$initialize(name, transformation)
    },
    calculate = function(seg_sites, trees, files, model) {
      ind <- get_population_individuals(model, private$population)
      calc_nucleotide_div(seg_sites, ind)
    }
  )
)

#' Summary Statistic: Nucleotide Diversity
#'
#' The summary statistic calculates the nucleotide diversity (\eqn{\pi})
#' per locus, which is the mean number of pairwise difference for
#' two individuals. It is a commonly used estimator for the scaled
#' mutation rate \eqn{\theta}.
#'
#' The nucleotide diversity was introduced by
#'
#'  Nei and Li (1979). "Mathematical Model for Studying Genetic Variation in
#'  Terms of Restriction Endonucleases". PNAS 76 (10): 5269-73.
#'  doi:10.1073/pnas.76.10.5269.
#'
#' @inheritParams sumstat_four_gamete
#' @return On simulation, this returns a vector with the value of pi for
#'   each locus.
#' @export
#' @template summary_statistics
#' @examples
#' model <- coal_model(5, 2) +
#'   feat_mutation(5) +
#'   sumstat_nucleotide_div()
#' stats <- simulate(model)
#' print(stats$pi)
sumstat_nucleotide_div <- function(name = "pi", population = 1,
                                   transformation = identity) {
  stat_pi_class$new(name, population, transformation)
}
