#' Identify clusters
#'
#' This is a special function used in the context of formula
#' used for Wilcoxon sum rank test for clustered data.
#' It identifies the cluster id of observations, and is used
#' on the right hand side of a formula.
#'
#' @param x A numeric variable of cluster id.
#'
#' @details THe function's only action is semantic, to mark
#' a variable as the cluster indicator. If not supplied,
#' will assume no cluster in the data.
#' @return x
#' @seealso cluswilcox.test.formula
#'
#' @examples
#' data(crd)
#' cluswilcox.test(z ~ cluster(id) + group(group), data = crd)
#' @export
cluster <- function(x) {x}



