#' @importFrom R6 R6Class
stat_jsfs_class <- R6Class("stat_jsfs", inherit = sumstat_class,
  private = list(
    populations = NULL,
    req_segsites = TRUE,
    per_locus = FALSE
  ),
  public = list(
    initialize = function(name, populations, per_locus, transformation) {
      assert_that(is.numeric(populations))
      private$populations <- populations

      assert_that(is.logical(per_locus))
      assert_that(length(per_locus) == 1)
      private$per_locus <- per_locus
      super$initialize(name, transformation)
    },
    calculate = function(seg_sites, trees, files, model) {
      ind_per_pop <- lapply(private$populations, get_population_individuals,
                            model = model)

      if (private$per_locus) {
        jsfs <- lapply(seg_sites, function(x) {
          calc_jsfs(list(x), ind_per_pop)
        })
      } else {
        jsfs <- calc_jsfs(seg_sites, ind_per_pop)
      }
      jsfs
    }
  )
)

#' Summary Statistic: Joint Site Frequency Spectrum
#'
#' The summary statistic calculates the joint site frequency
#' spectrum (JSFS) for multiple populations.
#'
#' @inheritParams sumstat_four_gamete
#' @param populations An integer vector containing the populations for which
#'        the JSFS is generated.
#' @param per_locus If \code{TRUE}, the JSFS is returned for each locus instead
#'   of globally. In this case, the result is a list, where each entry is the
#'   JSFS for the corresponding locus.
#' @return The JSFS, given as an array. The dimensions correspond to the
#'   populatons as given in the \code{populations} argument.
#' @template summary_statistics
#' @export
#' @examples
#' model <- coal_model(c(2, 3, 4), 2) +
#'   feat_mutation(5) +
#'   feat_migration(1, symmetric = TRUE) +
#'   sumstat_jsfs("jsfs_12", populations = c(1, 2)) +
#'   sumstat_jsfs("jsfs_123", populations = c(1, 2, 3))
#' stats <- simulate(model)
#' print(stats$jsfs_12)
#' print(stats$jsfs_123)
sumstat_jsfs <- function(name = "jsfs", populations = c(1, 2),
                         per_locus = FALSE, transformation = identity) {
  stat_jsfs_class$new(name, populations, per_locus, transformation)
}
