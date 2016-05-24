#'  Dataset with email messages and the folders where the users stored them
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 1702 instances, 1001 attributes and 53 labels
#' @source Klimt, B. and Yang, Y., "The Enron Corpus: A New Dataset for Email Classification Research", in Proc. ECML04, Pisa, Italy, pp. 217-226, 2004
#' @examples
#'\dontrun{
#' enron()  # Check and load the dataset
#' toBibtex(enron)
#' enron$measures
#' }
#' @export
enron <- function() check_n_load.mldr("enron")
