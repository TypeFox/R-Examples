#' @importFrom R6 R6Class
stat_tajimas_d_class <- R6Class("stat_tajimas_d", inherit = sumstat_class,
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
        n <- length(ind)
        if (n < 2) stop("At least two individuals are needed for Tajima's D")
        pi <- calc_nucleotide_div(seg_sites, ind)
        s <- sapply(seg_sites, ncol)

        tmp <- 1 / (1:(n - 1)) #nolint
        a1 <- sum(tmp)
        a2 <- sum(tmp ^ 2)
        b1 <- (n + 1) / (3 * (n - 1))
        b2 <- 2 * (n ^ 2 + n + 3) / (9 * n * (n - 1))
        c1 <- b1 - 1 / a1
        c2 <- b2 - (n + 2) / (a1 * n) + a2 / (a1 ^ 2)
        e1 <- c1 / a1
        e2 <- c2 / (a1 ^ 2 + a2)

        (pi - s / a1) / sqrt(e1 * s + e2 * s * (s - 1))
      }
   )
)


#' Summary Statistic: Tajima's D
#'
#' This statistic calculates Tajima's D from the simulation
#' results when added to a model. Tajima's D primarily measures
#' an deviation of singletons from the neutral expectation of
#' an equilibrium model. Negative values indicate an excess of
#' singletons, while positive values code a depletion of them.
#'
#' @inheritParams sumstat_four_gamete
#' @return On simulation, this returns a vector with the value of Tajima's D for
#'   each locus.
#' @references
#' Tajima, F. (1989). "Statistical method for testing the neutral mutation
#' hypothesis by DNA polymorphism.". Genetics 123 (3): 585-95.
#' @export
#' @template summary_statistics
#' @examples
#' # A neutral model that should yield values close to zero:
#' model <- coal_model(5, 2) +
#'   feat_mutation(20) +
#'   feat_recombination(10) +
#'   sumstat_tajimas_d("taji_d")
#' stats <- simulate(model)
#' print(stats$taji_d)
sumstat_tajimas_d <- function(name = "tajimas_d", population = "all",
                              transformation = identity) {
  stat_tajimas_d_class$new(name, population, transformation)
}
