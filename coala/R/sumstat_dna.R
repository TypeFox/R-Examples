#' @importFrom R6 R6Class
stat_dna_class <- R6Class("dna_stat", inherit = sumstat_class,
  private = list(),
  public = list(
    calculate = function(seg_sites, trees, dna, model) {
      if (requires_segsites(model)) {
        stop("Can not generate both seg. sites and DNA")
      }
      if (has_trios(model)) {
        stop("Cannot generate DNA for trio models")
      }
      lapply(dna, function(locus) {
        apply(locus, 2, function(x) attr(locus, "levels")[x])
      })
    }
  )
)

#' Summary Statistic: DNA
#'
#' This summary statistic returns the actual DNA sequences from
#' finite sites simulations. It can not be
#' calculated together with other summary statistics or when assuming
#' an infinite sites mutation model. No outgroup
#' is needed for it, and the outgroup sequences will also be
#' returned if present.
#'
#' @inheritParams sumstat_four_gamete
#' @return A list of sequences for each locus. Each entries is a
#'         character matrix decoding the sequences. Each row
#'         is an individual, and each column is a genetic position.
#' @template summary_statistics
#' @export
#' @examples
#' model <- coal_model(5, 1, 10) +
#'  feat_mutation(5, model = "GTR", gtr_rates = rep(1, 6)) +
#'  sumstat_dna()
#' \dontrun{simulate(model)$dna}
sumstat_dna <- function(name = "dna", transformation = identity) {
  stat_dna_class$new(name, transformation)
}
