#' @importFrom R6 R6Class
stat_segsites_class <- R6Class("stat_segsites", inherit = sumstat_class,
  private = list(req_segsites = TRUE),
  public = list(
    calculate = function(seg_sites, trees, files, model) seg_sites
  )
)


#' Summary Statistic: Segregating Sites
#'
#' This summary statistics generates a matrix of segregating sites.
#' This is useful for calculating summary statistics that \pkg{coala}
#' does not support..
#'
#' @inheritParams sumstat_four_gamete
#' @return A list of \code{\link[=create_segsites]{segsites}} objects.
#'         These can be treated as a matrix for
#'         most applications. Rows are individuals, columns are SNPs.
#' @export
#' @seealso For a description of the segregating sites class: \code{\link{create_segsites}}
#' @template summary_statistics
#' @examples
#' model <- coal_model(5, 1) +
#'   feat_mutation(5) +
#'   sumstat_seg_sites("segsites")
#' stats <- simulate(model)
#' print(stats$segsites)
sumstat_seg_sites <- function(name = "seg_sites", transformation = identity) {
  stat_segsites_class$new(name, transformation)
}


conv_for_trios <- function(seg_sites, model) {
  for (i in seq(along = seg_sites)) {
    locus_length <- get_locus_length(model, i, total = FALSE)
    if (length(locus_length) == 1) next

    total_length <- sum(locus_length)
    borders <- cumsum(locus_length[1:4] / total_length)

    pos <- get_positions(seg_sites[[i]])
    left <- pos < borders[1]
    middle <- pos >= borders[2] & pos < borders[3]
    right <- pos >= borders[4]

    pos[left] <- pos[left] * total_length / locus_length[1]
    pos[middle] <- (pos[middle] - borders[2]) * total_length / locus_length[3]
    pos[right] <- (pos[right] - borders[4]) * total_length / locus_length[5]

    trio_segsites <- seg_sites[[i]][ , left | middle | right]

    seg_sites[[i]] <- create_segsites(as.matrix(trio_segsites),
                                      pos[left | middle | right],
                                      c(rep(-1, sum(left)),
                                        rep(0, sum(middle)),
                                        rep(1, sum(right))),
                                      FALSE)

    assert_that(nrow(seg_sites[[i]]) > 0)
  }
  seg_sites
}
