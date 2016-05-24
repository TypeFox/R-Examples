#'  Dataset with protein profiles and their categories
#'
#' @description Multilabel dataset from the biology domain.
#' @format An mldr object with 2417 instances, 103 attributes and 14 labels
#' @source Elisseeff, A. and Weston, J., "A Kernel Method for Multi-Labelled Classification", Advances in Neural Information Processing Systems, Vol. 14, pp. 681--687, 2001
#' @examples
#'\dontrun{
#' yeast()  # Check and load the dataset
#' toBibtex(yeast)
#' yeast$measures
#' }
#' @export
yeast <- function() check_n_load.mldr('yeast')
