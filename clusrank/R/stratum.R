#' Identify stratums.
#'
#' This is a special function used in the context of formula
#' used for Wilcoxon sum rank test for clustered data.
#' It identifies the stratum id of observations, and is used
#' on the right hand side of a formula.
#'
#' @param x A numeric variable of stratum id.
#'
#' @details THe function's only action is semantic, to mark
#' a variable as the stratum indicator. If not supplied,
#' will assume no stratification in the data.
#' @seealso cluswilcox.test.formula
#'
#' @examples
#' data(crdStr)
#' cluswilcox.test(z ~ cluster(id) + group(group) + stratum(stratum), data = crdStr)
#' @export
stratum <- function(x) {x}