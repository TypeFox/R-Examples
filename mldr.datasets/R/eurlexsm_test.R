#'  List  with 10 folds of the test data from the EUR-Lex subject matters dataset
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 1935 instances, 5000 attributes and 201 labels
#' @source Mencia, E. L. and Furnkranz, J., "Efficient pairwise multilabel classification for large-scale problems in the legal domain", Machine Learning and Knowledge Discovery in Databases, pp. 50--65, 2008
#' @examples
#'\dontrun{
#' eurlexsm_test()  # Check and load the dataset
#' toBibtex(eurlexsm_test[[1]])
#' eurlexsm_test[[1]]$measures
#' }
#' @export
eurlexsm_test <- function() check_n_load.mldr('eurlexsm_test')
