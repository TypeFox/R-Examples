#'  Dataset generated from a subset of the Medline database
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 13929 instances, 1002 attributes and 23 labels
#' @source Joachims, Thorsten, "Text Categorization with Suport Vector Machines: Learning with Many Relevant Features", in Proc. 10th European Conference on Machine Learning, pp. 137-142, 1998
#' @examples
#'\dontrun{
#' ohsumed()  # Check and load the dataset
#' toBibtex(ohsumed)
#' ohsumed$measures
#' }
#' @export
ohsumed <- function() check_n_load.mldr('ohsumed')
