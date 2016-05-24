stat_trees_class <- R6Class("stat_trees", inherit = sumstat_class,
  private = list(req_files = FALSE, req_trees = TRUE),
  public = list(
    calculate = function(seg_sites, trees, files, model) {
      trees
    }
  )
)

#' Summary Statistic: Ancestral Trees
#'
#' This statistic returns ancestral tress in NEWICK format.
#'
#' @export
#' @inheritParams sumstat_four_gamete
#' @template summary_statistics
#' @examples
#' # Without recombination:
#' model <- coal_model(4, 2) + sumstat_trees()
#' stats <- simulate(model)
#' print(stats$trees)
#'
#' # With recombination:
#' model <- model + feat_recombination(5)
#' stats <- simulate(model)
#' print(stats$trees)
sumstat_trees <- function(name = "trees") {
  stat_trees_class$new(name, identity)
}


stat_sg_trees_class <- R6Class("stat_sg_trees", inherit = sumstat_class,
  private = list(req_trees = TRUE),
  public = list(
    calculate = function(seg_sites, trees, files, model) {
      llm <- get_locus_length_matrix(model)

      trio_trees <- generate_trio_trees(trees, llm)
      lapply(trio_trees, function(locus_trees) {
        files <- sapply(locus_trees, function(locus) {
          if (length(locus) > 0) {
            file <- tempfile("trio_tree")
            write(locus, file, sep = "\n")
          } else {
            file <- NULL
          }
          file
        })
        files[!sapply(files, is.null)]
      })
    }
  )
)

# Returns ancestral tress as files for seq-gen
sumstat_sg_trees <- function() {
  stat_sg_trees_class$new("trees", identity)
}
