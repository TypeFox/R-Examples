#' @importFrom R6 R6Class
stat_mcmf_class <- R6Class("stat_mcmf", inherit = sumstat_class,
  private = list(
    population = NULL,
    req_segsites = TRUE
  ),
  public = list(
    initialize = function(name, population, transformation) {
      assert_that(is.numeric(population))
      assert_that(length(population) == 1)
      private$population <- population
      super$initialize(name, transformation)
    },
    calculate = function(seg_sites, trees, files, model) {
      ploidy <- ifelse(is_unphased(model), get_ploidy(model), 1)
      calc_mcmf(seg_sites,
                get_population_individuals(model,
                                          private$population,
                                          haploids = (ploidy == 1)),
                has_trios(model),
                ploidy)
    }
  )
)

#' Summary Statistic: Most Common Mutation's Frequency
#'
#' This summary statistic calculates the observed frequency
#' of the mutational pattern that is observed most often in
#' the data.
#'
#' @inheritParams sumstat_four_gamete
#' @return A numeric vector containing MCMF for each locus.
#' @template summary_statistics
#' @examples
#' # Calculate MCMF for a panmictic population
#' model <- coal_model(10, 2) +
#'   feat_mutation(50) +
#'   sumstat_mcmf()
#' simulate(model)
#' @export
sumstat_mcmf  <- function(name = "mcmf", population = 1,
                          transformation = identity) {
  stat_mcmf_class$new(name, population, transformation)
}
