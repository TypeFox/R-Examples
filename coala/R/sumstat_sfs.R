#' @importFrom R6 R6Class
stat_sfs_class <- R6Class("stat_sfs", inherit = sumstat_class,
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
      individuals <- get_population_individuals(model, private$population)
      sfs <- as.vector(calc_jsfs(seg_sites, list(individuals)))
      sfs[c(-1, -length(sfs))]
    }
  )
)

#' Summary Statistic: Site Frequency Spectrum
#'
#' The Site Frequency Spectrum (SFS) counts how many
#' SNPs are in a sample according to their number of
#' derived alleles.
#'
#' @inheritParams sumstat_four_gamete
#' @export
#' @template summary_statistics
#' @examples
#' model <- coal_model(20, 500) +
#'   feat_mutation(2) +
#'   sumstat_sfs()
#' stats <- simulate(model)
#' barplot(stats$sfs)
sumstat_sfs <- function(name = "sfs", population = "all",
                        transformation = identity) {
  stat_sfs_class$new(name, population, transformation)
}
